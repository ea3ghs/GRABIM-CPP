TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -lnlopt

SOURCES += main.cpp \
    matchingnetwork.cpp

HEADERS += \
    matchingnetwork.h

QMAKE_CXXFLAGS += -DHAVE_STD -DHAVE_NAMESPACES -lopt -lnewmat -lnlopt -lm -std=c++11 -pthread -DNDEBUG
QMAKE_CXXFLAGS += -O2


#Manually: g++ -I/usr/local/include main.cpp matchingnetwork.cpp  -L/usr/local/lib -lnlopt -lm -o grabim
