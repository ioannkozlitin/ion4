#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QComboBox>
#include <QFileDialog>
#include <QMessageBox>
#include <QDebug>
#include <cmath>
#include <limits>

#define MIX_COMPONENTS_NUM 3

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

    for (int i = 0; i < MIX_COMPONENTS_NUM; i++)
    {
        on_addItemPushButton_clicked();
    }

    _ui->mixScrollArea->setFixedHeight(_ui->mixScrollAreaWidget->layout()->sizeHint().height());

    setFixedSize(0, 0);

    enableEditing(true);
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

    enableEditing(true);
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
    int newRow = layout->count() / 2;

    layout->addWidget(comboBox, newRow, 0);
    layout->addWidget(doubleSpinBox, newRow, 1);

    _ui->removeItemPushButton->setEnabled(newRow >= MIX_COMPONENTS_NUM);
}

void MainWindow::on_removeItemPushButton_clicked()
{
    QGridLayout *layout = (QGridLayout*)_ui->mixScrollAreaWidget->layout();
    int lastRow = layout->count() / 2 - 1;

    if (lastRow >= MIX_COMPONENTS_NUM)
    {
        disconnect(layout->itemAtPosition(lastRow, 0)->widget(), nullptr, nullptr, nullptr);

        layout->itemAtPosition(lastRow, 0)->widget()->hide();
        layout->itemAtPosition(lastRow, 1)->widget()->hide();

        layout->removeWidget(layout->itemAtPosition(lastRow, 0)->widget());
        layout->removeWidget(layout->itemAtPosition(lastRow, 1)->widget());
    }

    _ui->removeItemPushButton->setEnabled(lastRow > MIX_COMPONENTS_NUM);
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

        if (elementIndex > 0 && fabs(proportion) > std::numeric_limits<double>::min())
        {
            Z.push_back(elementIndex);
            x.push_back(proportion);
        }
    }

    if (Z.empty())
    {
        QMessageBox::warning(this, "Ошибка", "Не заданы компоненты смеси");
        return;
    }

    if (fabs(std::accumulate(x.begin(), x.end(), 0.0) - 1.0) > std::numeric_limits<double>::min())
    {
        QMessageBox::warning(this, "Ошибка", "Сумма долей компонентов смеси не равна единице");
        return;
    }

    double rCoeff = 0.6;
    double lgRhoMin = _ui->lgRhoMinSpinBox->value();
    double lgRhoMax = _ui->lgRhoMaxSpinBox->value();
    double lgRhoStep = _ui->lgRhoStepSpinBox->value();
    double lgTMin = _ui->lgTMinSpinBox->value();
    double lgTMax = _ui->lgTMaxSpinBox->value();
    double lgTStep = _ui->lgTStepSpinBox->value();

    if (lgRhoMax < lgRhoMin)
    {
        QMessageBox::warning(this, "Ошибка", "Неверно заданы lgRho");
        return;
    }

    if (lgTMax < lgTMin)
    {
        QMessageBox::warning(this, "Ошибка", "Неверно заданы lgT");
        return;
    }

    std::string filePath = _ui->filePathLineEdit->text().toStdString();

    if (filePath.empty())
    {
        QMessageBox::warning(this, "Ошибка", "Не выбран файл");
        return;
    }

    enableEditing(false);
    emit startCalc(Z, x, rCoeff, lgRhoMin, lgRhoMax, lgRhoStep, lgTMin, lgTMax, lgTStep, filePath);
}

void MainWindow::on_cancelPushButton_clicked()
{
    emit stopCalc();
}

void MainWindow::enableEditing(bool isEnabled)
{
    _ui->elementLabel->setEnabled(isEnabled);
    _ui->proportionLabel->setEnabled(isEnabled);
    _ui->mixScrollArea->setEnabled(isEnabled);
    _ui->addItemPushButton->setEnabled(isEnabled);

    int rowNum = ((QGridLayout*)_ui->mixScrollAreaWidget->layout())->count() / 2;
    _ui->removeItemPushButton->setEnabled(isEnabled && rowNum > MIX_COMPONENTS_NUM);

    _ui->minLabel->setEnabled(isEnabled);
    _ui->maxLabel->setEnabled(isEnabled);
    _ui->stepLabel->setEnabled(isEnabled);
    _ui->lgRhoLabel->setEnabled(isEnabled);
    _ui->lgRhoMinSpinBox->setEnabled(isEnabled);
    _ui->lgRhoMaxSpinBox->setEnabled(isEnabled);
    _ui->lgRhoStepSpinBox->setEnabled(isEnabled);
    _ui->lgTLabel->setEnabled(isEnabled);
    _ui->lgTMinSpinBox->setEnabled(isEnabled);
    _ui->lgTMaxSpinBox->setEnabled(isEnabled);
    _ui->lgTStepSpinBox->setEnabled(isEnabled);

    _ui->fileLabel->setEnabled(isEnabled);
    _ui->filePathLineEdit->setEnabled(isEnabled && !_ui->filePathLineEdit->text().isEmpty());
    _ui->selectFilePushButton->setEnabled(isEnabled);

    _ui->progressBar->setEnabled(!isEnabled);
    _ui->calculatePushButton->setEnabled(isEnabled);
    _ui->cancelPushButton->setEnabled(!isEnabled);
}
