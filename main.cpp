#include "dfwapcp.h"
#include <QtWidgets/QApplication>
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    DFWAPCP w;
    w.show();
    return a.exec();
}
