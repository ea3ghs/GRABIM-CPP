#include "ui.h"

using namespace std;


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    ui * GRABIM_interface =  new ui();
    GRABIM_interface->raise();
    GRABIM_interface->show();
    a.exec();
}



