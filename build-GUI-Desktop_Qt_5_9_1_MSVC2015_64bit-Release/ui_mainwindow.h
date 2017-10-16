/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.9.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QWidget>
#include "Output/qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QHBoxLayout *horizontalLayout;
    QWidget *widget;
    QTabWidget *tabWidget;
    QWidget *mainTab;
    QGroupBox *main_vcfGrp;
    QLineEdit *main_vcfDirTxt;
    QPushButton *main_vcfDirBtn;
    QCheckBox *main_vcfPassChk;
    QLineEdit *main_vcfMafTxt;
    QLabel *main_vcfMafLbl;
    QLabel *main_vcfMissingLbl;
    QLineEdit *main_vcfMissingTxt;
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout_3;
    QLineEdit *main_vcfChromFilterTxt;
    QLabel *main_vcfChrLbl;
    QPushButton *main_runBtn;
    QCustomPlot *main_plotBox;
    QGroupBox *main_sampleGrp;
    QPushButton *main_sampleDirBtn;
    QLineEdit *main_sampleDirTxt;
    QLineEdit *main_sampleDepthTxt;
    QLabel *main_sampleDepthLbl;
    QGroupBox *main_bedGrp;
    QPushButton *main_bedDirBtn;
    QLineEdit *main_bedDirTxt;
    QGroupBox *main_bedCollapseGrp;
    QHBoxLayout *horizontalLayout_2;
    QRadioButton *main_bedCollapseCodingBtn;
    QRadioButton *main_bedCollapseGeneBtn;
    QRadioButton *main_bedCollapseExonBtn;
    QGroupBox *main_testGrp;
    QRadioButton *main_testCommonBtn;
    QRadioButton *main_testRareCastBtn;
    QLineEdit *main_testBootTxt;
    QRadioButton *main_testRareCalphaBtn;
    QLabel *main_testBootLbl;
    QCheckBox *main_testBootChk;
    QWidget *simulatinTab;
    QGroupBox *sim_populationGrp;
    QGridLayout *gridLayout_7;
    QLineEdit *sim_populationSizeTxt;
    QLabel *sim_populationSizeLbl;
    QLineEdit *sim_populationPrevalenceTxt;
    QLabel *sim_pouulationPrevalenceLbl;
    QGroupBox *sim_sequencingGrp;
    QGridLayout *gridLayout_8;
    QLineEdit *sim_sequencingErrorRateTxt;
    QLabel *sim_sequencingErrorRateLbl;
    QLineEdit *sim_sequencingErrorSdTxt;
    QLabel *sim_sequencingErrorSdLbl;
    QGroupBox *sim_variantGrp;
    QLineEdit *sim_variantSizeTxt;
    QLabel *sim_variantMafLbl;
    QLineEdit *sim_variantMafTxt;
    QLabel *sim_variantSizeLbl;
    QTableWidget *sim_groupTbl;
    QPushButton *sim_runBtn;
    QGroupBox *sim_groupGrp;
    QGridLayout *gridLayout;
    QGroupBox *sim_groupDepthGrp;
    QLineEdit *sim_groupDepthMeanTxt;
    QLabel *sim_groupDepthMeanLbl;
    QLineEdit *sim_groupDepthSdTxt;
    QLabel *sim_groupDepthSdLbl;
    QGroupBox *sim_groupCohortGrp;
    QGridLayout *gridLayout_9;
    QRadioButton *sim_groupCaseBtn;
    QRadioButton *sim_groupControlBtn;
    QGroupBox *sim_groupSizeGrp;
    QLineEdit *sim_groupSizeMaxTxt;
    QLineEdit *sim_groupSizeMinTxt;
    QLabel *sim_groupDepthMeanLbl_4;
    QLabel *sim_groupDepthMeanLbl_5;
    QCustomPlot *sim_plotbox;
    QGroupBox *sim_testGrp;
    QRadioButton *sim_testCommonBtn;
    QRadioButton *sim_testRareCastBtn;
    QLineEdit *sim_testBootTxt;
    QRadioButton *sim_testRareCalphaBtn;
    QLabel *sim_testBootLbl;
    QCheckBox *sim_testBootChk;
    QPushButton *sim_groupAddBtn;
    QPushButton *sim_groupRemoveBtn;
    QLabel *sim_groupHighLowLbl;
    QLineEdit *sim_groupHighLowTxt;
    QGroupBox *sim_powerGrp;
    QLabel *sim_powerStepLbl;
    QLineEdit *sim_powerStepTxt;
    QLineEdit *sim_powerAlphaTxt;
    QLabel *sim_powerAlphaLbl;
    QLineEdit *sim_oddsRatioTxt;
    QLabel *sim_oddsRatioLbl;
    QWidget *tab;
    QLabel *label;
    QTextEdit *outputBox;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(919, 710);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        horizontalLayout = new QHBoxLayout(centralWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        widget = new QWidget(centralWidget);
        widget->setObjectName(QStringLiteral("widget"));
        tabWidget = new QTabWidget(widget);
        tabWidget->setObjectName(QStringLiteral("tabWidget"));
        tabWidget->setEnabled(true);
        tabWidget->setGeometry(QRect(0, 0, 885, 501));
        tabWidget->setMovable(true);
        tabWidget->setTabBarAutoHide(false);
        mainTab = new QWidget();
        mainTab->setObjectName(QStringLiteral("mainTab"));
        main_vcfGrp = new QGroupBox(mainTab);
        main_vcfGrp->setObjectName(QStringLiteral("main_vcfGrp"));
        main_vcfGrp->setGeometry(QRect(10, 0, 351, 131));
        main_vcfDirTxt = new QLineEdit(main_vcfGrp);
        main_vcfDirTxt->setObjectName(QStringLiteral("main_vcfDirTxt"));
        main_vcfDirTxt->setGeometry(QRect(10, 21, 251, 22));
        main_vcfDirTxt->setReadOnly(true);
        main_vcfDirBtn = new QPushButton(main_vcfGrp);
        main_vcfDirBtn->setObjectName(QStringLiteral("main_vcfDirBtn"));
        main_vcfDirBtn->setGeometry(QRect(280, 24, 20, 16));
        main_vcfDirBtn->setMaximumSize(QSize(20, 16777215));
        main_vcfPassChk = new QCheckBox(main_vcfGrp);
        main_vcfPassChk->setObjectName(QStringLiteral("main_vcfPassChk"));
        main_vcfPassChk->setGeometry(QRect(190, 60, 101, 16));
        main_vcfPassChk->setChecked(true);
        main_vcfMafTxt = new QLineEdit(main_vcfGrp);
        main_vcfMafTxt->setObjectName(QStringLiteral("main_vcfMafTxt"));
        main_vcfMafTxt->setGeometry(QRect(11, 61, 74, 22));
        main_vcfMafLbl = new QLabel(main_vcfGrp);
        main_vcfMafLbl->setObjectName(QStringLiteral("main_vcfMafLbl"));
        main_vcfMafLbl->setGeometry(QRect(89, 61, 91, 16));
        main_vcfMissingLbl = new QLabel(main_vcfGrp);
        main_vcfMissingLbl->setObjectName(QStringLiteral("main_vcfMissingLbl"));
        main_vcfMissingLbl->setGeometry(QRect(89, 90, 91, 16));
        main_vcfMissingTxt = new QLineEdit(main_vcfGrp);
        main_vcfMissingTxt->setObjectName(QStringLiteral("main_vcfMissingTxt"));
        main_vcfMissingTxt->setGeometry(QRect(11, 87, 74, 22));
        layoutWidget = new QWidget(main_vcfGrp);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(182, 84, 149, 24));
        horizontalLayout_3 = new QHBoxLayout(layoutWidget);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        main_vcfChromFilterTxt = new QLineEdit(layoutWidget);
        main_vcfChromFilterTxt->setObjectName(QStringLiteral("main_vcfChromFilterTxt"));

        horizontalLayout_3->addWidget(main_vcfChromFilterTxt);

        main_vcfChrLbl = new QLabel(layoutWidget);
        main_vcfChrLbl->setObjectName(QStringLiteral("main_vcfChrLbl"));

        horizontalLayout_3->addWidget(main_vcfChrLbl);

        main_runBtn = new QPushButton(mainTab);
        main_runBtn->setObjectName(QStringLiteral("main_runBtn"));
        main_runBtn->setGeometry(QRect(320, 400, 81, 21));
        main_plotBox = new QCustomPlot(mainTab);
        main_plotBox->setObjectName(QStringLiteral("main_plotBox"));
        main_plotBox->setGeometry(QRect(400, 10, 441, 331));
        main_sampleGrp = new QGroupBox(mainTab);
        main_sampleGrp->setObjectName(QStringLiteral("main_sampleGrp"));
        main_sampleGrp->setGeometry(QRect(10, 130, 351, 101));
        main_sampleDirBtn = new QPushButton(main_sampleGrp);
        main_sampleDirBtn->setObjectName(QStringLiteral("main_sampleDirBtn"));
        main_sampleDirBtn->setGeometry(QRect(280, 24, 20, 16));
        main_sampleDirBtn->setMaximumSize(QSize(20, 16777215));
        main_sampleDirTxt = new QLineEdit(main_sampleGrp);
        main_sampleDirTxt->setObjectName(QStringLiteral("main_sampleDirTxt"));
        main_sampleDirTxt->setGeometry(QRect(11, 21, 251, 22));
        main_sampleDirTxt->setReadOnly(true);
        main_sampleDepthTxt = new QLineEdit(main_sampleGrp);
        main_sampleDepthTxt->setObjectName(QStringLiteral("main_sampleDepthTxt"));
        main_sampleDepthTxt->setGeometry(QRect(10, 60, 51, 22));
        main_sampleDepthLbl = new QLabel(main_sampleGrp);
        main_sampleDepthLbl->setObjectName(QStringLiteral("main_sampleDepthLbl"));
        main_sampleDepthLbl->setGeometry(QRect(70, 60, 141, 22));
        main_bedGrp = new QGroupBox(mainTab);
        main_bedGrp->setObjectName(QStringLiteral("main_bedGrp"));
        main_bedGrp->setGeometry(QRect(10, 240, 341, 131));
        main_bedDirBtn = new QPushButton(main_bedGrp);
        main_bedDirBtn->setObjectName(QStringLiteral("main_bedDirBtn"));
        main_bedDirBtn->setGeometry(QRect(280, 30, 20, 16));
        main_bedDirBtn->setMaximumSize(QSize(20, 16777215));
        main_bedDirTxt = new QLineEdit(main_bedGrp);
        main_bedDirTxt->setObjectName(QStringLiteral("main_bedDirTxt"));
        main_bedDirTxt->setGeometry(QRect(10, 30, 251, 22));
        main_bedDirTxt->setReadOnly(true);
        main_bedCollapseGrp = new QGroupBox(main_bedGrp);
        main_bedCollapseGrp->setObjectName(QStringLiteral("main_bedCollapseGrp"));
        main_bedCollapseGrp->setGeometry(QRect(10, 60, 291, 61));
        horizontalLayout_2 = new QHBoxLayout(main_bedCollapseGrp);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        main_bedCollapseCodingBtn = new QRadioButton(main_bedCollapseGrp);
        main_bedCollapseCodingBtn->setObjectName(QStringLiteral("main_bedCollapseCodingBtn"));

        horizontalLayout_2->addWidget(main_bedCollapseCodingBtn);

        main_bedCollapseGeneBtn = new QRadioButton(main_bedCollapseGrp);
        main_bedCollapseGeneBtn->setObjectName(QStringLiteral("main_bedCollapseGeneBtn"));
        main_bedCollapseGeneBtn->setChecked(true);

        horizontalLayout_2->addWidget(main_bedCollapseGeneBtn);

        main_bedCollapseExonBtn = new QRadioButton(main_bedCollapseGrp);
        main_bedCollapseExonBtn->setObjectName(QStringLiteral("main_bedCollapseExonBtn"));

        horizontalLayout_2->addWidget(main_bedCollapseExonBtn);

        main_testGrp = new QGroupBox(mainTab);
        main_testGrp->setObjectName(QStringLiteral("main_testGrp"));
        main_testGrp->setGeometry(QRect(10, 370, 291, 87));
        main_testCommonBtn = new QRadioButton(main_testGrp);
        main_testCommonBtn->setObjectName(QStringLiteral("main_testCommonBtn"));
        main_testCommonBtn->setGeometry(QRect(10, 20, 91, 16));
        main_testCommonBtn->setChecked(true);
        main_testRareCastBtn = new QRadioButton(main_testGrp);
        main_testRareCastBtn->setObjectName(QStringLiteral("main_testRareCastBtn"));
        main_testRareCastBtn->setGeometry(QRect(10, 40, 111, 16));
        main_testBootTxt = new QLineEdit(main_testGrp);
        main_testBootTxt->setObjectName(QStringLiteral("main_testBootTxt"));
        main_testBootTxt->setEnabled(false);
        main_testBootTxt->setGeometry(QRect(140, 43, 74, 22));
        main_testRareCalphaBtn = new QRadioButton(main_testGrp);
        main_testRareCalphaBtn->setObjectName(QStringLiteral("main_testRareCalphaBtn"));
        main_testRareCalphaBtn->setGeometry(QRect(10, 60, 111, 16));
        main_testBootLbl = new QLabel(main_testGrp);
        main_testBootLbl->setObjectName(QStringLiteral("main_testBootLbl"));
        main_testBootLbl->setGeometry(QRect(218, 43, 71, 21));
        main_testBootChk = new QCheckBox(main_testGrp);
        main_testBootChk->setObjectName(QStringLiteral("main_testBootChk"));
        main_testBootChk->setGeometry(QRect(140, 23, 81, 16));
        tabWidget->addTab(mainTab, QString());
        simulatinTab = new QWidget();
        simulatinTab->setObjectName(QStringLiteral("simulatinTab"));
        sim_populationGrp = new QGroupBox(simulatinTab);
        sim_populationGrp->setObjectName(QStringLiteral("sim_populationGrp"));
        sim_populationGrp->setGeometry(QRect(10, 10, 141, 79));
        gridLayout_7 = new QGridLayout(sim_populationGrp);
        gridLayout_7->setSpacing(6);
        gridLayout_7->setContentsMargins(11, 11, 11, 11);
        gridLayout_7->setObjectName(QStringLiteral("gridLayout_7"));
        sim_populationSizeTxt = new QLineEdit(sim_populationGrp);
        sim_populationSizeTxt->setObjectName(QStringLiteral("sim_populationSizeTxt"));

        gridLayout_7->addWidget(sim_populationSizeTxt, 0, 0, 1, 1);

        sim_populationSizeLbl = new QLabel(sim_populationGrp);
        sim_populationSizeLbl->setObjectName(QStringLiteral("sim_populationSizeLbl"));

        gridLayout_7->addWidget(sim_populationSizeLbl, 0, 1, 1, 1);

        sim_populationPrevalenceTxt = new QLineEdit(sim_populationGrp);
        sim_populationPrevalenceTxt->setObjectName(QStringLiteral("sim_populationPrevalenceTxt"));

        gridLayout_7->addWidget(sim_populationPrevalenceTxt, 1, 0, 1, 1);

        sim_pouulationPrevalenceLbl = new QLabel(sim_populationGrp);
        sim_pouulationPrevalenceLbl->setObjectName(QStringLiteral("sim_pouulationPrevalenceLbl"));

        gridLayout_7->addWidget(sim_pouulationPrevalenceLbl, 1, 1, 1, 1);

        sim_sequencingGrp = new QGroupBox(simulatinTab);
        sim_sequencingGrp->setObjectName(QStringLiteral("sim_sequencingGrp"));
        sim_sequencingGrp->setGeometry(QRect(180, 10, 141, 79));
        gridLayout_8 = new QGridLayout(sim_sequencingGrp);
        gridLayout_8->setSpacing(6);
        gridLayout_8->setContentsMargins(11, 11, 11, 11);
        gridLayout_8->setObjectName(QStringLiteral("gridLayout_8"));
        sim_sequencingErrorRateTxt = new QLineEdit(sim_sequencingGrp);
        sim_sequencingErrorRateTxt->setObjectName(QStringLiteral("sim_sequencingErrorRateTxt"));

        gridLayout_8->addWidget(sim_sequencingErrorRateTxt, 0, 0, 1, 1);

        sim_sequencingErrorRateLbl = new QLabel(sim_sequencingGrp);
        sim_sequencingErrorRateLbl->setObjectName(QStringLiteral("sim_sequencingErrorRateLbl"));

        gridLayout_8->addWidget(sim_sequencingErrorRateLbl, 0, 1, 1, 1);

        sim_sequencingErrorSdTxt = new QLineEdit(sim_sequencingGrp);
        sim_sequencingErrorSdTxt->setObjectName(QStringLiteral("sim_sequencingErrorSdTxt"));

        gridLayout_8->addWidget(sim_sequencingErrorSdTxt, 1, 0, 1, 1);

        sim_sequencingErrorSdLbl = new QLabel(sim_sequencingGrp);
        sim_sequencingErrorSdLbl->setObjectName(QStringLiteral("sim_sequencingErrorSdLbl"));

        gridLayout_8->addWidget(sim_sequencingErrorSdLbl, 1, 1, 1, 1);

        sim_variantGrp = new QGroupBox(simulatinTab);
        sim_variantGrp->setObjectName(QStringLiteral("sim_variantGrp"));
        sim_variantGrp->setGeometry(QRect(10, 380, 151, 81));
        sim_variantSizeTxt = new QLineEdit(sim_variantGrp);
        sim_variantSizeTxt->setObjectName(QStringLiteral("sim_variantSizeTxt"));
        sim_variantSizeTxt->setGeometry(QRect(9, 22, 51, 22));
        sim_variantMafLbl = new QLabel(sim_variantGrp);
        sim_variantMafLbl->setObjectName(QStringLiteral("sim_variantMafLbl"));
        sim_variantMafLbl->setGeometry(QRect(70, 50, 31, 16));
        sim_variantMafTxt = new QLineEdit(sim_variantGrp);
        sim_variantMafTxt->setObjectName(QStringLiteral("sim_variantMafTxt"));
        sim_variantMafTxt->setGeometry(QRect(9, 48, 51, 22));
        sim_variantSizeLbl = new QLabel(sim_variantGrp);
        sim_variantSizeLbl->setObjectName(QStringLiteral("sim_variantSizeLbl"));
        sim_variantSizeLbl->setGeometry(QRect(70, 22, 71, 16));
        sim_groupTbl = new QTableWidget(simulatinTab);
        if (sim_groupTbl->columnCount() < 4)
            sim_groupTbl->setColumnCount(4);
        QTableWidgetItem *__qtablewidgetitem = new QTableWidgetItem();
        sim_groupTbl->setHorizontalHeaderItem(0, __qtablewidgetitem);
        QTableWidgetItem *__qtablewidgetitem1 = new QTableWidgetItem();
        sim_groupTbl->setHorizontalHeaderItem(1, __qtablewidgetitem1);
        QTableWidgetItem *__qtablewidgetitem2 = new QTableWidgetItem();
        sim_groupTbl->setHorizontalHeaderItem(2, __qtablewidgetitem2);
        QTableWidgetItem *__qtablewidgetitem3 = new QTableWidgetItem();
        sim_groupTbl->setHorizontalHeaderItem(3, __qtablewidgetitem3);
        sim_groupTbl->setObjectName(QStringLiteral("sim_groupTbl"));
        sim_groupTbl->setEnabled(true);
        sim_groupTbl->setGeometry(QRect(10, 240, 321, 141));
        sim_groupTbl->horizontalHeader()->setStretchLastSection(true);
        sim_groupTbl->verticalHeader()->setStretchLastSection(false);
        sim_runBtn = new QPushButton(simulatinTab);
        sim_runBtn->setObjectName(QStringLiteral("sim_runBtn"));
        sim_runBtn->setGeometry(QRect(730, 420, 80, 21));
        sim_runBtn->setAutoDefault(false);
        sim_runBtn->setFlat(false);
        sim_groupGrp = new QGroupBox(simulatinTab);
        sim_groupGrp->setObjectName(QStringLiteral("sim_groupGrp"));
        sim_groupGrp->setGeometry(QRect(10, 90, 311, 106));
        gridLayout = new QGridLayout(sim_groupGrp);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        sim_groupDepthGrp = new QGroupBox(sim_groupGrp);
        sim_groupDepthGrp->setObjectName(QStringLiteral("sim_groupDepthGrp"));
        sim_groupDepthMeanTxt = new QLineEdit(sim_groupDepthGrp);
        sim_groupDepthMeanTxt->setObjectName(QStringLiteral("sim_groupDepthMeanTxt"));
        sim_groupDepthMeanTxt->setGeometry(QRect(10, 23, 41, 20));
        sim_groupDepthMeanLbl = new QLabel(sim_groupDepthGrp);
        sim_groupDepthMeanLbl->setObjectName(QStringLiteral("sim_groupDepthMeanLbl"));
        sim_groupDepthMeanLbl->setGeometry(QRect(60, 23, 26, 16));
        sim_groupDepthSdTxt = new QLineEdit(sim_groupDepthGrp);
        sim_groupDepthSdTxt->setObjectName(QStringLiteral("sim_groupDepthSdTxt"));
        sim_groupDepthSdTxt->setGeometry(QRect(10, 49, 41, 20));
        sim_groupDepthSdLbl = new QLabel(sim_groupDepthGrp);
        sim_groupDepthSdLbl->setObjectName(QStringLiteral("sim_groupDepthSdLbl"));
        sim_groupDepthSdLbl->setGeometry(QRect(60, 49, 16, 16));

        gridLayout->addWidget(sim_groupDepthGrp, 0, 0, 1, 1);

        sim_groupCohortGrp = new QGroupBox(sim_groupGrp);
        sim_groupCohortGrp->setObjectName(QStringLiteral("sim_groupCohortGrp"));
        gridLayout_9 = new QGridLayout(sim_groupCohortGrp);
        gridLayout_9->setSpacing(6);
        gridLayout_9->setContentsMargins(11, 11, 11, 11);
        gridLayout_9->setObjectName(QStringLiteral("gridLayout_9"));
        sim_groupCaseBtn = new QRadioButton(sim_groupCohortGrp);
        sim_groupCaseBtn->setObjectName(QStringLiteral("sim_groupCaseBtn"));
        sim_groupCaseBtn->setChecked(true);

        gridLayout_9->addWidget(sim_groupCaseBtn, 0, 0, 1, 1);

        sim_groupControlBtn = new QRadioButton(sim_groupCohortGrp);
        sim_groupControlBtn->setObjectName(QStringLiteral("sim_groupControlBtn"));

        gridLayout_9->addWidget(sim_groupControlBtn, 1, 0, 1, 1);


        gridLayout->addWidget(sim_groupCohortGrp, 0, 1, 1, 1);

        sim_groupSizeGrp = new QGroupBox(sim_groupGrp);
        sim_groupSizeGrp->setObjectName(QStringLiteral("sim_groupSizeGrp"));
        sim_groupSizeMaxTxt = new QLineEdit(sim_groupSizeGrp);
        sim_groupSizeMaxTxt->setObjectName(QStringLiteral("sim_groupSizeMaxTxt"));
        sim_groupSizeMaxTxt->setGeometry(QRect(10, 49, 31, 20));
        sim_groupSizeMinTxt = new QLineEdit(sim_groupSizeGrp);
        sim_groupSizeMinTxt->setObjectName(QStringLiteral("sim_groupSizeMinTxt"));
        sim_groupSizeMinTxt->setGeometry(QRect(10, 23, 31, 20));
        sim_groupDepthMeanLbl_4 = new QLabel(sim_groupSizeGrp);
        sim_groupDepthMeanLbl_4->setObjectName(QStringLiteral("sim_groupDepthMeanLbl_4"));
        sim_groupDepthMeanLbl_4->setGeometry(QRect(50, 23, 16, 16));
        sim_groupDepthMeanLbl_5 = new QLabel(sim_groupSizeGrp);
        sim_groupDepthMeanLbl_5->setObjectName(QStringLiteral("sim_groupDepthMeanLbl_5"));
        sim_groupDepthMeanLbl_5->setGeometry(QRect(50, 49, 20, 16));

        gridLayout->addWidget(sim_groupSizeGrp, 0, 2, 1, 1);

        sim_groupCohortGrp->raise();
        sim_groupDepthGrp->raise();
        sim_groupSizeGrp->raise();
        sim_plotbox = new QCustomPlot(simulatinTab);
        sim_plotbox->setObjectName(QStringLiteral("sim_plotbox"));
        sim_plotbox->setGeometry(QRect(430, 50, 391, 241));
        sim_testGrp = new QGroupBox(simulatinTab);
        sim_testGrp->setObjectName(QStringLiteral("sim_testGrp"));
        sim_testGrp->setGeometry(QRect(170, 380, 231, 81));
        sim_testCommonBtn = new QRadioButton(sim_testGrp);
        sim_testCommonBtn->setObjectName(QStringLiteral("sim_testCommonBtn"));
        sim_testCommonBtn->setGeometry(QRect(10, 20, 70, 16));
        sim_testCommonBtn->setChecked(true);
        sim_testRareCastBtn = new QRadioButton(sim_testGrp);
        sim_testRareCastBtn->setObjectName(QStringLiteral("sim_testRareCastBtn"));
        sim_testRareCastBtn->setGeometry(QRect(10, 40, 91, 16));
        sim_testBootTxt = new QLineEdit(sim_testGrp);
        sim_testBootTxt->setObjectName(QStringLiteral("sim_testBootTxt"));
        sim_testBootTxt->setEnabled(false);
        sim_testBootTxt->setGeometry(QRect(110, 50, 51, 22));
        sim_testRareCalphaBtn = new QRadioButton(sim_testGrp);
        sim_testRareCalphaBtn->setObjectName(QStringLiteral("sim_testRareCalphaBtn"));
        sim_testRareCalphaBtn->setGeometry(QRect(10, 60, 91, 16));
        sim_testBootLbl = new QLabel(sim_testGrp);
        sim_testBootLbl->setObjectName(QStringLiteral("sim_testBootLbl"));
        sim_testBootLbl->setGeometry(QRect(170, 50, 61, 16));
        sim_testBootChk = new QCheckBox(sim_testGrp);
        sim_testBootChk->setObjectName(QStringLiteral("sim_testBootChk"));
        sim_testBootChk->setGeometry(QRect(110, 23, 91, 16));
        sim_groupAddBtn = new QPushButton(simulatinTab);
        sim_groupAddBtn->setObjectName(QStringLiteral("sim_groupAddBtn"));
        sim_groupAddBtn->setGeometry(QRect(10, 210, 80, 21));
        sim_groupRemoveBtn = new QPushButton(simulatinTab);
        sim_groupRemoveBtn->setObjectName(QStringLiteral("sim_groupRemoveBtn"));
        sim_groupRemoveBtn->setGeometry(QRect(100, 210, 91, 21));
        sim_groupHighLowLbl = new QLabel(simulatinTab);
        sim_groupHighLowLbl->setObjectName(QStringLiteral("sim_groupHighLowLbl"));
        sim_groupHighLowLbl->setGeometry(QRect(280, 210, 91, 22));
        sim_groupHighLowTxt = new QLineEdit(simulatinTab);
        sim_groupHighLowTxt->setObjectName(QStringLiteral("sim_groupHighLowTxt"));
        sim_groupHighLowTxt->setGeometry(QRect(200, 210, 74, 22));
        sim_powerGrp = new QGroupBox(simulatinTab);
        sim_powerGrp->setObjectName(QStringLiteral("sim_powerGrp"));
        sim_powerGrp->setGeometry(QRect(410, 380, 261, 81));
        sim_powerStepLbl = new QLabel(sim_powerGrp);
        sim_powerStepLbl->setObjectName(QStringLiteral("sim_powerStepLbl"));
        sim_powerStepLbl->setGeometry(QRect(90, 22, 41, 16));
        sim_powerStepTxt = new QLineEdit(sim_powerGrp);
        sim_powerStepTxt->setObjectName(QStringLiteral("sim_powerStepTxt"));
        sim_powerStepTxt->setGeometry(QRect(9, 22, 74, 22));
        sim_powerAlphaTxt = new QLineEdit(sim_powerGrp);
        sim_powerAlphaTxt->setObjectName(QStringLiteral("sim_powerAlphaTxt"));
        sim_powerAlphaTxt->setGeometry(QRect(9, 48, 74, 22));
        sim_powerAlphaLbl = new QLabel(sim_powerGrp);
        sim_powerAlphaLbl->setObjectName(QStringLiteral("sim_powerAlphaLbl"));
        sim_powerAlphaLbl->setGeometry(QRect(90, 48, 31, 16));
        sim_oddsRatioTxt = new QLineEdit(sim_powerGrp);
        sim_oddsRatioTxt->setObjectName(QStringLiteral("sim_oddsRatioTxt"));
        sim_oddsRatioTxt->setGeometry(QRect(130, 40, 51, 22));
        sim_oddsRatioLbl = new QLabel(sim_powerGrp);
        sim_oddsRatioLbl->setObjectName(QStringLiteral("sim_oddsRatioLbl"));
        sim_oddsRatioLbl->setGeometry(QRect(190, 40, 71, 22));
        tabWidget->addTab(simulatinTab, QString());
        tab = new QWidget();
        tab->setObjectName(QStringLiteral("tab"));
        label = new QLabel(tab);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(80, 30, 121, 81));
        tabWidget->addTab(tab, QString());
        outputBox = new QTextEdit(widget);
        outputBox->setObjectName(QStringLiteral("outputBox"));
        outputBox->setGeometry(QRect(0, 510, 881, 151));

        horizontalLayout->addWidget(widget);

        MainWindow->setCentralWidget(centralWidget);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        MainWindow->setStatusBar(statusBar);

        retranslateUi(MainWindow);

        tabWidget->setCurrentIndex(0);
        sim_runBtn->setDefault(false);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "vikNGS GUI interface", Q_NULLPTR));
        main_vcfGrp->setTitle(QApplication::translate("MainWindow", "Multisample VCF File", Q_NULLPTR));
        main_vcfDirTxt->setText(QApplication::translate("MainWindow", "C:/Users/Scott/Desktop/RVS-master/example/example_1000snps.vcf", Q_NULLPTR));
        main_vcfDirBtn->setText(QApplication::translate("MainWindow", "...", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        main_vcfPassChk->setToolTip(QApplication::translate("MainWindow", "QUAL column in VCF must contain PASS", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        main_vcfPassChk->setText(QApplication::translate("MainWindow", "Must PASS", Q_NULLPTR));
        main_vcfMafTxt->setText(QApplication::translate("MainWindow", "0.05", Q_NULLPTR));
        main_vcfMafTxt->setPlaceholderText(QApplication::translate("MainWindow", "0.05", Q_NULLPTR));
        main_vcfMafLbl->setText(QApplication::translate("MainWindow", "MAF Cutoff", Q_NULLPTR));
        main_vcfMissingLbl->setText(QApplication::translate("MainWindow", "Missing Threshold", Q_NULLPTR));
        main_vcfMissingTxt->setText(QApplication::translate("MainWindow", "0.2", Q_NULLPTR));
        main_vcfChromFilterTxt->setText(QString());
        main_vcfChrLbl->setText(QApplication::translate("MainWindow", "Chromosome filter", Q_NULLPTR));
        main_runBtn->setText(QApplication::translate("MainWindow", "RUN", Q_NULLPTR));
        main_sampleGrp->setTitle(QApplication::translate("MainWindow", "Sample Data File", Q_NULLPTR));
        main_sampleDirBtn->setText(QApplication::translate("MainWindow", "...", Q_NULLPTR));
        main_sampleDirTxt->setText(QString());
        main_sampleDepthTxt->setText(QApplication::translate("MainWindow", "30", Q_NULLPTR));
        main_sampleDepthLbl->setText(QApplication::translate("MainWindow", "Read Depth High/Low Cutoff", Q_NULLPTR));
        main_bedGrp->setTitle(QApplication::translate("MainWindow", "BED file (optional)", Q_NULLPTR));
        main_bedDirBtn->setText(QApplication::translate("MainWindow", "...", Q_NULLPTR));
        main_bedCollapseGrp->setTitle(QApplication::translate("MainWindow", "Collapse By", Q_NULLPTR));
        main_bedCollapseCodingBtn->setText(QApplication::translate("MainWindow", "Coding only", Q_NULLPTR));
        main_bedCollapseGeneBtn->setText(QApplication::translate("MainWindow", "Entire gene", Q_NULLPTR));
        main_bedCollapseExonBtn->setText(QApplication::translate("MainWindow", "Exon only", Q_NULLPTR));
        main_testGrp->setTitle(QApplication::translate("MainWindow", "Test", Q_NULLPTR));
        main_testCommonBtn->setText(QApplication::translate("MainWindow", "Common", Q_NULLPTR));
        main_testRareCastBtn->setText(QApplication::translate("MainWindow", "Rare (CAST)", Q_NULLPTR));
        main_testBootTxt->setText(QApplication::translate("MainWindow", "100", Q_NULLPTR));
        main_testRareCalphaBtn->setText(QApplication::translate("MainWindow", "Rare (c-alpha)", Q_NULLPTR));
        main_testBootLbl->setText(QApplication::translate("MainWindow", "# iterations", Q_NULLPTR));
        main_testBootChk->setText(QApplication::translate("MainWindow", "Bootstrap", Q_NULLPTR));
        tabWidget->setTabText(tabWidget->indexOf(mainTab), QApplication::translate("MainWindow", "Association Test", Q_NULLPTR));
        sim_populationGrp->setTitle(QApplication::translate("MainWindow", "Population Simulation", Q_NULLPTR));
        sim_populationSizeTxt->setText(QApplication::translate("MainWindow", "1000", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_populationSizeLbl->setToolTip(QApplication::translate("MainWindow", "Size of population to simulate", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_populationSizeLbl->setText(QApplication::translate("MainWindow", "Size", Q_NULLPTR));
        sim_populationPrevalenceTxt->setText(QApplication::translate("MainWindow", "0.2", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_pouulationPrevalenceLbl->setToolTip(QApplication::translate("MainWindow", "Prevalence of the disease in the simulated population", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_pouulationPrevalenceLbl->setText(QApplication::translate("MainWindow", "Prevalence", Q_NULLPTR));
        sim_sequencingGrp->setTitle(QApplication::translate("MainWindow", "Sequencing Parameters", Q_NULLPTR));
        sim_sequencingErrorRateTxt->setText(QApplication::translate("MainWindow", "0.01", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_sequencingErrorRateLbl->setToolTip(QApplication::translate("MainWindow", "Base calling error rate of sequencer", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_sequencingErrorRateLbl->setText(QApplication::translate("MainWindow", "Error Rate", Q_NULLPTR));
        sim_sequencingErrorSdTxt->setText(QApplication::translate("MainWindow", "0.025", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_sequencingErrorSdLbl->setToolTip(QApplication::translate("MainWindow", "Standard deviation of base calling error rate of sequencer", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_sequencingErrorSdLbl->setText(QApplication::translate("MainWindow", "Error SD", Q_NULLPTR));
        sim_variantGrp->setTitle(QApplication::translate("MainWindow", "Variant Parameters", Q_NULLPTR));
        sim_variantSizeTxt->setText(QApplication::translate("MainWindow", "100", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_variantMafLbl->setToolTip(QApplication::translate("MainWindow", "Minor allele frequency", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_variantMafLbl->setText(QApplication::translate("MainWindow", "MAF", Q_NULLPTR));
        sim_variantMafTxt->setText(QApplication::translate("MainWindow", "0.1", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_variantSizeLbl->setToolTip(QApplication::translate("MainWindow", "Number of variants to test", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_variantSizeLbl->setText(QApplication::translate("MainWindow", "# of Variants", Q_NULLPTR));
        QTableWidgetItem *___qtablewidgetitem = sim_groupTbl->horizontalHeaderItem(0);
        ___qtablewidgetitem->setText(QApplication::translate("MainWindow", "# Individuals", Q_NULLPTR));
        QTableWidgetItem *___qtablewidgetitem1 = sim_groupTbl->horizontalHeaderItem(1);
        ___qtablewidgetitem1->setText(QApplication::translate("MainWindow", "Cohort", Q_NULLPTR));
        QTableWidgetItem *___qtablewidgetitem2 = sim_groupTbl->horizontalHeaderItem(2);
        ___qtablewidgetitem2->setText(QApplication::translate("MainWindow", "Mean Depth", Q_NULLPTR));
        QTableWidgetItem *___qtablewidgetitem3 = sim_groupTbl->horizontalHeaderItem(3);
        ___qtablewidgetitem3->setText(QApplication::translate("MainWindow", "Depth SD", Q_NULLPTR));
        sim_runBtn->setText(QApplication::translate("MainWindow", "RUN", Q_NULLPTR));
        sim_groupGrp->setTitle(QApplication::translate("MainWindow", "Group", Q_NULLPTR));
        sim_groupDepthGrp->setTitle(QApplication::translate("MainWindow", "Sequencing Depth", Q_NULLPTR));
        sim_groupDepthMeanTxt->setText(QApplication::translate("MainWindow", "30", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_groupDepthMeanLbl->setToolTip(QApplication::translate("MainWindow", "Mean read depth of sequencing group", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_groupDepthMeanLbl->setText(QApplication::translate("MainWindow", "Mean", Q_NULLPTR));
        sim_groupDepthSdTxt->setText(QApplication::translate("MainWindow", "5", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_groupDepthSdLbl->setToolTip(QApplication::translate("MainWindow", "Read depth standard deviation of sequencing group", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_groupDepthSdLbl->setText(QApplication::translate("MainWindow", "SD", Q_NULLPTR));
        sim_groupCohortGrp->setTitle(QApplication::translate("MainWindow", "Cohort", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_groupCaseBtn->setToolTip(QApplication::translate("MainWindow", "Group contains affected individuals", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_groupCaseBtn->setText(QApplication::translate("MainWindow", "Case", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_groupControlBtn->setToolTip(QApplication::translate("MainWindow", "Group contains unaffected individuals", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_groupControlBtn->setText(QApplication::translate("MainWindow", "Control", Q_NULLPTR));
        sim_groupSizeGrp->setTitle(QApplication::translate("MainWindow", "Size", Q_NULLPTR));
        sim_groupSizeMaxTxt->setText(QApplication::translate("MainWindow", "40", Q_NULLPTR));
        sim_groupSizeMinTxt->setText(QApplication::translate("MainWindow", "20", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_groupDepthMeanLbl_4->setToolTip(QApplication::translate("MainWindow", "Mean read depth of sequencing group", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_groupDepthMeanLbl_4->setText(QApplication::translate("MainWindow", "Min", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_groupDepthMeanLbl_5->setToolTip(QApplication::translate("MainWindow", "Mean read depth of sequencing group", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_groupDepthMeanLbl_5->setText(QApplication::translate("MainWindow", "Max", Q_NULLPTR));
        sim_testGrp->setTitle(QApplication::translate("MainWindow", "Test", Q_NULLPTR));
        sim_testCommonBtn->setText(QApplication::translate("MainWindow", "Common", Q_NULLPTR));
        sim_testRareCastBtn->setText(QApplication::translate("MainWindow", "Rare (CAST)", Q_NULLPTR));
        sim_testBootTxt->setText(QApplication::translate("MainWindow", "100", Q_NULLPTR));
        sim_testRareCalphaBtn->setText(QApplication::translate("MainWindow", "Rare (c-alpha)", Q_NULLPTR));
        sim_testBootLbl->setText(QApplication::translate("MainWindow", "# iterations", Q_NULLPTR));
        sim_testBootChk->setText(QApplication::translate("MainWindow", "Bootstrap", Q_NULLPTR));
        sim_groupAddBtn->setText(QApplication::translate("MainWindow", "Add Group", Q_NULLPTR));
        sim_groupRemoveBtn->setText(QApplication::translate("MainWindow", "Remove Selected", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_groupHighLowLbl->setToolTip(QApplication::translate("MainWindow", "Any mean read depth less than this will be considered a low read group", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_groupHighLowLbl->setText(QApplication::translate("MainWindow", "High/low cutoff", Q_NULLPTR));
        sim_groupHighLowTxt->setText(QApplication::translate("MainWindow", "30", Q_NULLPTR));
        sim_powerGrp->setTitle(QApplication::translate("MainWindow", "Power", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_powerStepLbl->setToolTip(QApplication::translate("MainWindow", "Mean read depth of sequencing group", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_powerStepLbl->setText(QApplication::translate("MainWindow", "Steps", Q_NULLPTR));
        sim_powerStepTxt->setText(QApplication::translate("MainWindow", "20", Q_NULLPTR));
        sim_powerAlphaTxt->setText(QApplication::translate("MainWindow", "0.05", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_powerAlphaLbl->setToolTip(QApplication::translate("MainWindow", "Odds ratio of disease", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_powerAlphaLbl->setText(QApplication::translate("MainWindow", "Alpha", Q_NULLPTR));
        sim_oddsRatioTxt->setText(QApplication::translate("MainWindow", "1.5", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        sim_oddsRatioLbl->setToolTip(QApplication::translate("MainWindow", "Odds ratio of disease", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        sim_oddsRatioLbl->setText(QApplication::translate("MainWindow", "Odds Ratio", Q_NULLPTR));
        tabWidget->setTabText(tabWidget->indexOf(simulatinTab), QApplication::translate("MainWindow", "Power Simulation (case/control)", Q_NULLPTR));
        label->setText(QApplication::translate("MainWindow", "Coming Soon!", Q_NULLPTR));
        tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("MainWindow", "Power Simulation (quantitative)", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
