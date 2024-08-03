/*****************************************************************************
 * Alpine Terrain Renderer
 * Copyright (C) 2024 Adam Celarek
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "TrackModel.h"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#else
#include <QFileDialog>
#endif

#include <gl_engine/Context.h>

TrackModel::TrackModel(QObject* parent)
    : QObject { parent }
{

    auto& c = gl_engine::Context::instance();
    c.setup_tracks(&m_manager);
}

#ifdef __EMSCRIPTEN__
// clang-format off
EM_ASYNC_JS(void, alpine_app_open_file_picker_and_mount, (), {
    const file_reader = new FileReader();
    const fileReadPromise = new Promise((resolve, reject) => {
        file_reader.onload = (event) => {
            try {
                const uint8Arr = new Uint8Array(event.target.result);
                const stream = FS.open("/tmp/track_upload/" + file_reader.filename.replace(/[^a-zA-Z0-9_\\-()]/g, '_'), "w+");
                FS.write(stream, uint8Arr, 0, uint8Arr.length, 0);
                FS.close(stream);
                resolve();
            } catch (error) {
                reject(error);
            }
        };
    });
    globalThis["open_file"] = function(e)
    {
        file_reader.filename = e.target.files[0].name;
        file_reader.mime_type = e.target.files[0].type;
        file_reader.readAsArrayBuffer(e.target.files[0]);
    };
    
    var file_selector = document.createElement('input');
    file_selector.setAttribute('type', 'file');
    file_selector.setAttribute('onchange', 'globalThis["open_file"](event)');
    file_selector.setAttribute('accept', '.gpx');
    file_selector.click();
    await fileReadPromise;
});
// clang-format on
#endif

void TrackModel::upload_track()
{
    auto fileContentReady = [this](const QString& /*fileName*/, const QByteArray& fileContent) {
        (void)fileContent;
        QXmlStreamReader xmlReader(fileContent);

        std::unique_ptr<nucleus::gpx::Gpx> gpx = nucleus::gpx::parse(xmlReader);
        if (gpx != nullptr) {
            m_manager.add_or_replace(0, *gpx);

            // if (0 < gpx->track.size() && 0 < gpx->track[0].size()) {
            //     auto track_start = gpx->track[0][0];
            //     emit position_set_by_user(track_start.latitude, track_start.longitude);
            // }
        } else {
            qDebug("Coud not parse GPX file!");
        }
    };

#ifdef __EMSCRIPTEN__
    QDir dir("/tmp/track_upload/");
    dir.mkpath("/tmp/track_upload/");
    alpine_app_open_file_picker_and_mount();
    QStringList fileList = dir.entryList(QDir::Files);

    if (fileList.isEmpty()) {
        qDebug() << "No files in /tmp/track_upload/";
        return;
    }
    QString file_name = fileList.first();
    QFile file(dir.absoluteFilePath(file_name));

    if (!file.open(QIODevice::ReadOnly)) {
        qDebug() << "Failed to open the file:" << file_name;
        return;
    }

    QByteArray file_data = file.readAll();
    file.close();

    if (!QFile::remove(dir.absoluteFilePath(file_name))) {
        qDebug() << "Failed to delete the file:" << file_name;
    }
    fileContentReady(file_name, file_data);
#else
    QFileDialog::getOpenFileContent("GPX (*.gpx *.xml)", fileContentReady);
#endif
}
