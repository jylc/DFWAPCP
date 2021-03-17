#include "imgcanvas.h"
#include <qevent.h>
imgcanvas::imgcanvas(QWidget* parent)
	:QLabel(parent), m_nMsgCount(0), m_pParent(parent)
{
    setAlignment(Qt::AlignCenter);
	m_nMovePointX = 0;
	m_nMovePointY = 0;
	m_Timer.setInterval(0);
	connect(&m_Timer, SIGNAL(QTimer::timeout()), this, SLOT(msgCountAdd()));
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
        qDebug() << ev->globalPos() - dragPosition;
        m_nMovePointX = (ev->globalPos() - dragPosition).x();
        m_nMovePointY = (ev->globalPos() - dragPosition).y();
        if (m_nMovePointX < 0 || m_nMovePointY < 0 || m_nMovePointY > m_pParent->rect().bottom() - this->rect().width() || m_nMovePointX > m_pParent->rect().right() - this->rect().width())
        {
            //��ֹ�Ƴ������͵ײ�
            if (m_nMovePointY > m_pParent->rect().bottom() - this->rect().width() || m_nMovePointY < 0)
            {
                move((ev->globalPos() - dragPosition).x(), this->geometry().y());
            }
            //��ֹ�Ƴ������Ҳ�
            if (m_nMovePointX > m_pParent->rect().right() - this->rect().width() || m_nMovePointX < 0)
            {
                move(this->geometry().x(), (ev->globalPos() - dragPosition).y());
            }
            //��ֹ�Ƴ����Ͻ�
            if (m_nMovePointX < 0 && m_nMovePointY < 0)
            {
                move(0, 0);
            }
            //��ֹ�Ƴ���С��
            if (m_nMovePointY > m_pParent->rect().bottom() - this->rect().width() && m_nMovePointX > m_pParent->rect().right() - this->rect().width())
            {
                move(m_pParent->rect().right() - this->rect().width(), m_pParent->rect().bottom() - this->rect().height());
            }
            //��ֹ�Ƴ����Ͻ�
            if (m_nMovePointY < 0 && m_nMovePointX > m_pParent->rect().right() - this->rect().width())
            {
                move(m_pParent->rect().right() - this->rect().width(), 0);
            }
            //��ֹ�Ƴ����Ͻ�
            if (m_nMovePointX < 0 && m_nMovePointY > m_pParent->rect().bottom() - this->rect().height())
            {
                move(0, m_pParent->rect().bottom() - this->rect().height());
            }

        }
        else
        {
            move(ev->globalPos() - dragPosition);
        }
    }
}
void imgcanvas::mouseReleaseEvent(QMouseEvent* ev)
{

}