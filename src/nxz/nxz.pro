QT += core
QT -= gui

TARGET = nxz
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

INCLUDEPATH += ../../../vcglib ../../../vcglib/eigenlib

SOURCES += main.cpp \
    nxzdecoder.cpp \
    nxzencoder.cpp
HEADERS += \
    nxzdecoder.h \
    nxzencoder.h

DISTFILES += \
    plan.md

