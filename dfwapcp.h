#pragma once

#include <QtWidgets/QWidget>
#include "ui_dfwapcp.h"

class DFWAPCP : public QWidget
{
    Q_OBJECT

public:
    DFWAPCP(QWidget *parent = Q_NULLPTR);

private:
    Ui::PanoWindow1 ui;
};
