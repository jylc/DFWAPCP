#include "imgcanvas.h"
#include <qevent.h>
imgcanvas::imgcanvas(QWidget* parent)
	:QLabel(parent), m_nMsgCount(0), m_pParent(parent)
{
    setAlignment(Qt::AlignCenter);
	m_nMovePointX = 0;
	m_nMovePointY = 0;
	m_Timer.setInterval(0);
	connect(&m_Timer, SIGNAL(timeout()), this, SLOT(msgCountAdd()));
}

void imgcanvas::msgCountAdd()
{
    update();
}
void imgcanvas::mousePressEvent(QMouseEvent* ev)
{
	if (ev->button() == Qt::LeftButton)
	{
		dragPosition = ev->globalPos() - frameGeometry().topLeft();
		ev->accept();
	}
}
void imgcanvas::mouseMoveEvent(QMouseEvent* ev)
{
    if (ev->buttons() & Qt::LeftButton)
    {
		move(ev->globalPos() - dragPosition);
    }
}
void imgcanvas::mouseReleaseEvent(QMouseEvent* ev)
{

}