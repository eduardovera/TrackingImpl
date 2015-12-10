TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += include

SOURCES += main.cpp

HEADERS += \
    settings.h \
    utils.h \
    Point.h \
    ransac.h

LIBS += -lX11 -lpthread
