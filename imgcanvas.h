#pragma once
#include <qlabel.h>
#include <qtimer.h>
#include <qdebug.h>
class imgcanvas :
    public QLabel
{
    Q_OBJECT
public:
    explicit imgcanvas(QWidget* parent = 0);

signals:
protected:
    void mousePressEvent(QMouseEvent* ev);
    void mouseMoveEvent(QMouseEvent* ev);
    void mouseReleaseEvent(QMouseEvent* ev);

public slots:
    void msgCountAdd();

private:
    QPoint dragPosition;
    QWidget* m_pParent;
    int m_nMovePointX;
    int m_nMovePointY;
    int m_nMsgCount;
    QTimer m_Timer;
};

