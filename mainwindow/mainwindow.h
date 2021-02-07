#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui { class MainWindow; }

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
    void setCalcProgress(int percentage);
    void calcFinished(bool isSuccess);

private slots:
    void on_addItemPushButton_clicked();
    void on_removeItemPushButton_clicked();
    void on_selectFilePushButton_clicked();
    void on_calculatePushButton_clicked();
    void on_cancelPushButton_clicked();

signals:
    void startCalc(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep, const std::string &filePath);
    void stopCalc();

private:
    Ui::MainWindow *_ui;

    QStringList _elements;
    QStringList _sortedElements;
};

#endif // MAINWINDOW_H
