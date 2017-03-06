QT += core
QT -= gui

TARGET = nxz
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++11

INCLUDEPATH += ../../../vcglib ../../../vcglib/eigenlib

SOURCES += main.cpp \
    nxzdecoder.cpp \
    nxzencoder.cpp \
    tunstall.cpp \
    bitstream.cpp \
    ../../../vcglib/wrap/ply/plylib.cpp
HEADERS += \
    nxzdecoder.h \
    nxzencoder.h \
    point.h \
    zpoint.h \
    cstream.h \
    tunstall.h \
    bitstream.h \
    ../../../vcglib/wrap/ply/plylib.h \
    nxz.h

DISTFILES += \
    plan.md

