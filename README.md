# AlpineMaps.org Renderer
This is the software behind [alpinemaps.org](https://alpinemaps.org).

A developer version (trunk) is released [here](https://alpinemapsorg.github.io/renderer/), including APKs for android. Be aware that it can break at any time!

[If looking at the issues, best to filter out projects!](https://github.com/AlpineMapsOrg/renderer/issues?q=is%3Aissue%20state%3Aopen%20no%3Aproject)

We are in discord, talk to us!
https://discord.gg/p8T9XzVwRa

# Cloning and building
`git clone git@github.com:AlpineMapsOrg/renderer.git`

After that it should be a normal cmake project. That is, you run cmake to generate a project or build file and then run your favourite tool. All dependencies should be pulled automatically while you run CMake. 
We use Qt Creator (with mingw on Windows), which is the only tested setup atm and makes setup of Android and WebAssembly builds reasonably easy. If you have questions, please go to Discord.

## Dependencies
* Qt 6.8.0, or greater
* g++ 12+, clang or msvc
* OpenGL
* Qt Positioning and Charts modules
* Some other dependencies will be pulled automatically during building.

## Building the native version
* just run cmake and build

## Building the android version
* We are usually building with Qt Creator, because it works relatively out of the box. However, it should also work on the command line or other IDEs if you set it up correctly.
* You need a Java JDK before you can do anything else. Not all Java versions work, and the error messages might be surprising (or non-existant). I'm running with Java 19, and I can compile for old devices. Iirc a newer version of Java caused issues. [Android documents the required Java version](https://developer.android.com/build/jdks), but as said, for me Java 19 works as well. It might change in the future.
* Once you have Java, go to Qt Creator Preferences -> Devices -> Android. There click "Set Up SDK" to automatically download and install an Android SDK.
* Finally, you might need to click on SDK Manager to install a fitting SDK Platform (take the newest, it also works for older devices), and ndk (newest as well).
* Then Google the internet to find out how to enable the developer mode on Android.
* On linux, you'll have to setup some udev rules. Run `Android/SDK/platform-tools/adb devices` and you should get instructions.
* If there are problems, check out the [documentation from Qt](https://doc.qt.io/qt-6/android-getting-started.html) 
* Finally, you are welcome to ask in discord if something is not working! 

## Building the WebAssembly version:
* [The Qt documentation is quite good on how to get it to run](https://doc-snapshots.qt.io/qt6-dev/wasm.html#installing-emscripten).
* Be aware that only specific versions of emscripten work for specific versions of Qt, and the error messages are not helpfull.
* [More info on building and getting Hotreload to work](https://github.com/AlpineMapsOrg/documentation/blob/main/WebAssembly_local_build.md)

# Code style
* class names are CamelCase, method, function and variable names are snake_case.
* class attributes have an m_ prefix and are usually private, struct attributes don't and are usually public.
* use `void set_attribute(int value)` and `int attribute() const` for setters and getters (that is, avoid the get_). Use [the Qt recommendations](https://wiki.qt.io/API_Design_Principles#Naming_Boolean_Getters,_Setters,_and_Properties) for naming boolean getters.
* structs are usually small, simple, and have no or only few methods. they never have inheritance.
* files are CamelCase if the content is a CamelCase class. otherwise they are snake_case, and have a snake_case namespace with stuff.
* the folder/structure.h is reflected in namespace folder::structure{ .. }
* indent with space only, indent 4 spaces
* ideally, use the clang-format file provided with the project
  (in case you use Qt Creator, go to Preferences -> C++ -> Code Style: Formatting mode: Full, Format while typing, Format edited code on file save, don't override formatting)
* follow the [Qt recommendations](https://wiki.qt.io/API_Design_Principles) and the [c++ core guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines) everywhere else.

# Developer workflow
* Fork this repository.
* Enable github pages from actions (Repository Settings -> Pages -> Source -> GitHub Actions)
* Work in branches or your main.
* Make pull requests from your main branch.
* Github Actions will run the unit tests and create packages for the browser and Android and deploy them to your_username.github.io/your_clone_name/.
* Make sure that the unit tests run through.
* We will also look at the browser version during the pull request.
* Ideally you'll also setup the signing keys for Android packages ([instructions](https://github.com/AlpineMapsOrg/renderer/blob/main/creating_apk_keys.md)).
