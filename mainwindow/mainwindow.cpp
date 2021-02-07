#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QComboBox>
#include <QFileDialog>
#include <QMessageBox>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , _ui(new Ui::MainWindow)
{
    _ui->setupUi(this);

    _elements = QStringList({""  , "H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne", "Na", "Mg", "Al", "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca", "Sc", "Ti", "V" , "Cr", "Mn",
                             "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
                             "Te", "I" , "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir",
                             "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U" , "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"});

    _sortedElements = _elements;
    qSort(_sortedElements);

    for (int i = 0; i < 3; i++)
    {
        on_addItemPushButton_clicked();
    }
}

MainWindow::~MainWindow()
{
    delete _ui;
}

void MainWindow::setCalcProgress(int percentage)
{
    _ui->progressBar->setValue(percentage);
}

void MainWindow::calcFinished(bool isSuccess)
{
    _ui->progressBar->setValue(0);

    if (isSuccess)
    {
        QMessageBox::information(this, "Сообщение", "Результат вычислений сохранен в файл");
    }
}

void MainWindow::on_addItemPushButton_clicked()
{
    QComboBox *comboBox = new QComboBox;
    comboBox->addItems(_sortedElements);
    comboBox->setStyleSheet("combobox-popup: 0");
    comboBox->setMaxVisibleItems(10);

    QDoubleSpinBox *doubleSpinBox = new QDoubleSpinBox;
    doubleSpinBox->setRange(0.0, 1.0);
    doubleSpinBox->setDecimals(6);
    doubleSpinBox->setSingleStep(1e-6);
    doubleSpinBox->setEnabled(false);

    connect(comboBox, &QComboBox::currentTextChanged, doubleSpinBox, [doubleSpinBox] (const QString &text) { doubleSpinBox->setEnabled(text != ""); });

    QGridLayout *layout = (QGridLayout*)_ui->mixScrollAreaWidget->layout();
    int row = layout->count() / 2;

    layout->addWidget(comboBox, row, 0);
    layout->addWidget(doubleSpinBox, row, 1);
}

void MainWindow::on_removeItemPushButton_clicked()
{
    QGridLayout *layout = (QGridLayout*)_ui->mixScrollAreaWidget->layout();
    int lastRow = layout->count() / 2 - 1;

    if (lastRow > 2)
    {
        disconnect(layout->itemAtPosition(lastRow, 0)->widget(), nullptr, nullptr, nullptr);

        layout->itemAtPosition(lastRow, 0)->widget()->hide();
        layout->itemAtPosition(lastRow, 1)->widget()->hide();

        layout->removeWidget(layout->itemAtPosition(lastRow, 0)->widget());
        layout->removeWidget(layout->itemAtPosition(lastRow, 1)->widget());
    }
}

void MainWindow::on_selectFilePushButton_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this, "Сохранение файла", "", "*.m");

    if (fileName.size() != 0)
    {
        if (!fileName.endsWith(".m"))
        {
            fileName += ".m";
        }

        _ui->filePathLineEdit->setEnabled(true);
        _ui->filePathLineEdit->setText(fileName);
    }
}

void MainWindow::on_calculatePushButton_clicked()
{
    std::vector<unsigned int> Z;
    std::vector<double> x;

    QGridLayout *layout = (QGridLayout*)_ui->mixScrollAreaWidget->layout();

    for (int i = 0; i < layout->count() / 2; i++)
    {
        QComboBox *comboBox = (QComboBox*)layout->itemAtPosition(i, 0)->widget();
        int elementIndex = _elements.indexOf(comboBox->currentText());

        QDoubleSpinBox *doubleSpinBox = (QDoubleSpinBox*)layout->itemAtPosition(i, 1)->widget();
        double proportion = doubleSpinBox->value();

        if (elementIndex > 0)
        {
            Z.push_back(elementIndex);
            x.push_back(proportion);
        }        
    }

    double rCoeff = 0.6;
    double lgRhoMin = _ui->lgRhoMinSpinBox->value();
    double lgRhoMax = _ui->lgRhoMaxSpinBox->value();
    double lgRhoStep = _ui->lgRhoStepSpinBox->value();
    double lgTMin = _ui->lgTMinSpinBox->value();
    double lgTMax = _ui->lgTMaxSpinBox->value();
    double lgTStep = _ui->lgTStepSpinBox->value();
    std::string filePath = _ui->filePathLineEdit->text().toStdString();

    emit startCalc(Z, x, rCoeff, lgRhoMin, lgRhoMax, lgRhoStep, lgTMin, lgTMax, lgTStep, filePath);
}

void MainWindow::on_cancelPushButton_clicked()
{
    emit stopCalc();
}
