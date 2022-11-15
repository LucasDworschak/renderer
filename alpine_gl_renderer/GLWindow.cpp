/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include <array>
#include <iostream>

#include <QDebug>
#include <QImage>
#include <QMoveEvent>
#include <QOpenGLBuffer>
#include <QOpenGLContext>
#include <QOpenGLDebugLogger>
#include <QOpenGLExtraFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>
#include <QOpenGLVertexArrayObject>
#include <QPropertyAnimation>
#include <QRandomGenerator>
#include <QSequentialAnimationGroup>
#include <QTimer>

#include <glm/glm.hpp>

#include "GLWindow.h"
#include "alpine_gl_renderer/GLDebugPainter.h"
#include "alpine_gl_renderer/GLShaderManager.h"
#include "alpine_gl_renderer/GLTileManager.h"
#include "nucleus/TileScheduler.h"

GLWindow::GLWindow()
    : m_camera({ 1822577.0, 6141664.0 - 500, 171.28 + 500 }, { 1822577.0, 6141664.0, 171.28 }) // should point right at the stephansdom
{
    QTimer::singleShot(0, [this]() { this->update(); });
}

GLWindow::~GLWindow()
{
    makeCurrent();
}
constexpr auto depth_internal_format = GL_DEPTH_COMPONENT32F;
constexpr auto depth_type = GL_FLOAT;
//constexpr auto depth_internal_format = GL_DEPTH_COMPONENT16;
//constexpr auto depth_type = GL_UNSIGNED_SHORT;
//constexpr auto depth_internal_format = GL_DEPTH_COMPONENT32;
//constexpr auto depth_type = GL_INT;

void GLWindow::initializeGL()
{
    QOpenGLDebugLogger* logger = new QOpenGLDebugLogger(this);
    logger->initialize();
    connect(logger, &QOpenGLDebugLogger::messageLogged, [](const auto& message) {
        qDebug() << message;
    });
    logger->disableMessages(QList<GLuint>({ 131185 }));
    logger->startLogging(QOpenGLDebugLogger::SynchronousLogging);
    const auto c = QOpenGLContext::currentContext();
    QOpenGLFunctions* f = QOpenGLContext::currentContext()->extraFunctions();

    m_tile_manager = std::make_unique<GLTileManager>();
    m_debug_painter = std::make_unique<GLDebugPainter>();
    m_shader_manager = std::make_unique<GLShaderManager>();

    m_tile_manager->setAttributeLocations(m_shader_manager->tileAttributeLocations());
    m_tile_manager->setUniformLocations(m_shader_manager->tileUniformLocations());

    m_debug_painter->setAttributeLocations(m_shader_manager->debugAttributeLocations());
    m_debug_painter->setUniformLocations(m_shader_manager->debugUniformLocations());

    {
        f->glGenTextures(1, &m_frame_buffer_colour);
        f->glBindTexture(GL_TEXTURE_2D, m_frame_buffer_colour);
        {
            const auto level = 0;
            const auto internalFormat = GL_RGBA;
            const auto border = 0;
            const auto format = GL_RGBA;
            const auto type = GL_UNSIGNED_BYTE;
            const auto data = nullptr;
            f->glTexImage2D(GL_TEXTURE_2D, level, internalFormat, 4, 4, border, format, type, data);
        }
    }
    {
        f->glGenTextures(1, &m_frame_buffer_depth);
        f->glBindTexture(GL_TEXTURE_2D, m_frame_buffer_depth);
        {
            const auto level = 0;
            const auto internalFormat = depth_internal_format;
            const auto border = 0;
            const auto format = GL_DEPTH_COMPONENT;
            const auto type = depth_type;
            const auto data = nullptr;
            f->glTexImage2D(GL_TEXTURE_2D, level, internalFormat, 4, 4, border, format, type, data);
        }
    }
    f->glBindTexture(GL_TEXTURE_2D, 0);
    {
        f->glGenFramebuffers(1, &m_frame_buffer);
        f->glBindFramebuffer(GL_FRAMEBUFFER, m_frame_buffer);
        f->glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_frame_buffer_colour, 0);
        f->glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, m_frame_buffer_depth, 0);
    }
//    std::cout << "init complete" << std::endl;
}

void GLWindow::resizeGL(int w, int h)
{
    if (w == 0 || h == 0)
        return;
    const qreal retinaScale = devicePixelRatio();
    const int width = int(retinaScale * w);
    const int height = int(retinaScale * h);

    QOpenGLFunctions* f = QOpenGLContext::currentContext()->functions();
    {
        f->glBindTexture(GL_TEXTURE_2D, m_frame_buffer_colour);
        const auto level = 0;
        const auto internalFormat = GL_RGBA;
        const auto border = 0;
        const auto format = GL_RGBA;
        const auto type = GL_UNSIGNED_BYTE;
        const auto data = nullptr;
        f->glTexImage2D(GL_TEXTURE_2D, level, internalFormat, width, height, border, format, type, data);
    }
    {
        f->glBindTexture(GL_TEXTURE_2D, m_frame_buffer_depth);
        const auto level = 0;
        const auto internalFormat = depth_internal_format;
        const auto border = 0;
        const auto format = GL_DEPTH_COMPONENT;
        const auto type = depth_type;
        const auto data = nullptr;
        f->glTexImage2D(GL_TEXTURE_2D, level, internalFormat, width, height, border, format, type, data);
    }

    f->glViewport(0, 0, width, height);
    emit viewport_changed({ w, h });
//    std::cout << "resizeGL" << std::endl;
}

void GLWindow::paintGL()
{
//    std::cout << "paintGL start" << std::endl;
    m_frame_start = std::chrono::time_point_cast<ClockResolution>(Clock::now());

    QOpenGLExtraFunctions* f = QOpenGLContext::currentContext()->extraFunctions();
    f->glBindFramebuffer(GL_FRAMEBUFFER, m_frame_buffer);
    f->glClearColor(1.0, 0.0, 0.5, 1);

    f->glClearDepthf(1.00);
    f->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    f->glEnable(GL_DEPTH_TEST);
    f->glDepthFunc(GL_LESS);
//    f->glDepthFunc(GL_GREATER);
//    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

    m_shader_manager->bindTileShader();

    const auto world_view_projection_matrix = m_camera.localViewProjectionMatrix({});
    m_tile_manager->draw(m_shader_manager->tileShader(), world_view_projection_matrix);

    {
        m_shader_manager->bindDebugShader();
        m_debug_painter->activate(m_shader_manager->debugShader(), world_view_projection_matrix);
        const auto position = m_camera.position();
        const auto direction_tl = m_camera.ray_direction({ -1, 1 });
        const auto direction_tr = m_camera.ray_direction({ 1, 1 });
        std::vector<glm::vec3> debug_cam_lines = { position + direction_tl * 10000.0,
            position,
            position + direction_tr * 10000.0 };
        m_debug_painter->drawLineStrip(debug_cam_lines);
    }
    m_shader_manager->release();

//        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

//        f->glDepthFunc(GL_LESS);
//        f->glDisable(GL_DEPTH_TEST);

    const qreal retinaScale = devicePixelRatio();
    const int width = int(retinaScale * this->width());
    const int height = int(retinaScale * this->height());
    f->glBindFramebuffer(GL_READ_FRAMEBUFFER, m_frame_buffer);
    f->glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    f->glBlitFramebuffer(0, 0, width, height, 0, 0, width, height, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    f->glBindFramebuffer(GL_FRAMEBUFFER, 0);
    m_frame_end = std::chrono::time_point_cast<ClockResolution>(Clock::now());

//    QTimer::singleShot(500, [this]() { this->update(); });
//    std::cout << "paintGL end" << std::endl;
}

void GLWindow::paintOverGL()
{
    const auto frame_duration = (m_frame_end - m_frame_start);
    const auto frame_duration_float = double(frame_duration.count()) / 1000.;
    const auto frame_duration_text = QString("Last frame: %1ms, draw indicator: ")
                                         .arg(QString::asprintf("%04.1f", frame_duration_float));

    const auto scheduler_stats = QString("Scheduler: %1 tiles in transit, %2 waiting height tiles, %3 waiting ortho tiles, %4 tiles on gpu")
                                     .arg(m_tile_scheduler->numberOfTilesInTransit())
                                     .arg(m_tile_scheduler->numberOfWaitingHeightTiles())
                                     .arg(m_tile_scheduler->numberOfWaitingOrthoTiles())
                                     .arg(m_tile_scheduler->gpuTiles().size());

    const auto random_u32 = QRandomGenerator::global()->generate();

    QPainter painter(this);
    painter.setFont(QFont("Helvetica", 12));
    painter.setPen(Qt::white);
    QRect text_bb = painter.boundingRect(10, 20, 1, 15, Qt::TextSingleLine, frame_duration_text);
    painter.drawText(10, 20, frame_duration_text);
    painter.drawText(10, 40, scheduler_stats);
    painter.setBrush(QBrush(QColor(random_u32)));
    painter.drawRect(int(text_bb.right()) + 5, 8, 12, 12);
}

void GLWindow::mouseMoveEvent(QMouseEvent* e)
{
    // send depth information only on mouse press to be more efficient (?)
    emit mouse_moved(e);
}

void GLWindow::keyPressEvent(QKeyEvent* e)
{
    if (e->key() == Qt::Key::Key_T) {
        m_tile_scheduler->setEnabled(!m_tile_scheduler->enabled());
        qDebug("setting tile scheduler enabled = %d", int(m_tile_scheduler->enabled()));
    }
    emit key_pressed(e);
}

void GLWindow::update_camera(const camera::Definition& new_definition)
{
    m_camera = new_definition;
    update();
}

void GLWindow::setTileScheduler(TileScheduler* new_tile_scheduler)
{
    m_tile_scheduler = new_tile_scheduler;
}

void GLWindow::mousePressEvent(QMouseEvent* ev)
{
    float distance = 500;   // todo read and compute from depth buffer
    emit mouse_pressed(ev, distance);
}

GLTileManager* GLWindow::gpuTileManager() const
{
    return m_tile_manager.get();
}
