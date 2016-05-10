#include "ui.h"

ui::ui()
{
      setFixedSize(400, 500);
      centralWidget =  new QWidget();
      QHBoxLayout *Impedancelayout = new QHBoxLayout();
      QHBoxLayout * ButtonsLayout = new QHBoxLayout();
      QGridLayout * TopoLayout = new QGridLayout();
      QHBoxLayout * LocalOptLayout = new QHBoxLayout();

      //Create source load impedance buttons
      SourceFileButton = new QPushButton("Source impedance");
      LoadFileButton = new QPushButton("Load impedance");
      connect(SourceFileButton, SIGNAL(clicked()), this, SLOT(SourceImpedance_clicked()));
      connect(LoadFileButton, SIGNAL(clicked()), this, SLOT(LoadImpedance_clicked()));


      QVBoxLayout * SourceLayout = new QVBoxLayout();
      QVBoxLayout * LoadLayout = new QVBoxLayout();

      QHBoxLayout * ConstantZSLayout = new QHBoxLayout();
      QHBoxLayout * ConstantZLLayout = new QHBoxLayout();

      FixedZSLineedit = new QLineEdit("50");
      FixedZLLineedit = new QLineEdit("50");
      ZSOhmLabel = new QLabel("Ohm");
      ZLOhmLabel = new QLabel("Ohm");

      ZSOhmLabel->setVisible(false);
      ZLOhmLabel->setVisible(false);

      FixedZSLineedit->setVisible(false);
      SourceFileButton->setVisible(true);

      FixedZLLineedit->setVisible(false);
      LoadFileButton->setVisible(true);

      ConstantZSLayout->addWidget(FixedZSLineedit);
      ConstantZSLayout->addWidget(ZSOhmLabel);

      ConstantZLLayout->addWidget(FixedZLLineedit);
      ConstantZLLayout->addWidget(ZLOhmLabel);

      SourceLayout->addWidget(SourceFileButton);
      SourceLayout->addLayout(ConstantZSLayout);

      LoadLayout->addWidget(LoadFileButton);
      LoadLayout->addLayout(ConstantZLLayout);

      //Checkboxes for selecting a constant impedance vs freq profile
      FixedZSCheckbox = new QCheckBox("Constant ZS");
      FixedZLCheckbox = new QCheckBox("Constant ZL");
      connect(FixedZSCheckbox, SIGNAL(clicked()), this, SLOT(FixedZSCheckbox_clicked()));
      connect(FixedZLCheckbox, SIGNAL(clicked()), this, SLOT(FixedZLCheckbox_clicked()));

      SourceLayout->addWidget(FixedZSCheckbox);
      LoadLayout->addWidget(FixedZLCheckbox);

      // Image
       QSize sz;
       QString s1 = "./../GRABIM-CPP/bitmaps/MatchingNetwork.svg";
       imgWidget = new QSvgWidget(s1);
       sz = imgWidget->size();
       imgWidget->setFixedSize(.2*sz);


      Impedancelayout->addLayout(SourceLayout);
      Impedancelayout->addWidget(imgWidget);
      Impedancelayout->addLayout(LoadLayout);




       // Matching band
       QGroupBox *groupBox = new QGroupBox(tr("Matching band"));
       minFLabel = new QLabel("Min:");
       minFEdit = new QLineEdit("1");
       minFUnitsCombo =  new QComboBox();
       maxFLabel = new QLabel("Max:");
       maxFEdit = new QLineEdit("2");
       maxFUnitsCombo =  new QComboBox();


       // Fill combos
       minFUnitsCombo->addItem("kHz");
       minFUnitsCombo->addItem("MHz");
       minFUnitsCombo->addItem("GHz");
       minFUnitsCombo->setCurrentIndex(2);

       maxFUnitsCombo->addItem("kHz");
       maxFUnitsCombo->addItem("MHz");
       maxFUnitsCombo->addItem("GHz");
       maxFUnitsCombo->setCurrentIndex(2);

       QHBoxLayout *vbox = new QHBoxLayout;
       vbox->addWidget(minFLabel);
       vbox->addWidget(minFEdit);
       vbox->addWidget(minFUnitsCombo);

       vbox->addWidget(maxFLabel);
       vbox->addWidget(maxFEdit);
       vbox->addWidget(maxFUnitsCombo);
       groupBox->setLayout(vbox);

       //Topology selection
       ArbitraryTopology = new QRadioButton("User-defined topology");
       ArbitraryTopology->setChecked(false);
       SearchBestTopology = new QRadioButton("Search best topology");
       SearchBestTopology->setChecked(true);

       QHBoxLayout *helpbox = new QHBoxLayout;
       CodeLabel = new QLabel("0: Series inductance\n1: Series capacitor\n2: Paralell inductance\n3: Parallel capacitor\n4: Transmission line\n5: Open stub\n6: Short circuited stub");
       CodeLabel->setWordWrap(true);
       CodeLabel->setVisible(false);
       helpbox->addWidget(CodeLabel);

       ArbitraryTopologyLineEdit =  new QLineEdit();

       ArbitraryTopologyLineEdit->setEnabled(false);
       TopoLayout->addWidget(ArbitraryTopology, 0, 0);
       TopoLayout->addWidget(ArbitraryTopologyLineEdit, 0, 1);
       TopoLayout->addWidget(SearchBestTopology, 1, 0);


       connect(ArbitraryTopology, SIGNAL(clicked()), this, SLOT(ArbitraryTopology_clicked()));
       connect(SearchBestTopology, SIGNAL(clicked()), this, SLOT(SearchBestTopology_clicked()));


       //Local optimizer settings
       LocalOptLabel = new QLabel("Local optimiser");
       LocalOptCombo =  new QComboBox();
       LocalOptCombo->addItem("Nelder-Mead");
       LocalOptCombo->addItem("PRAXIS");
       LocalOptCombo->addItem("Suplex");
       LocalOptCombo->addItem("COBYLA");
       LocalOptCombo->addItem("BOBYQA");
       LocalOptCombo->addItem("Augmented Lagrangian");
       LocalOptCombo->addItem("ESCH");
       LocalOptCombo->addItem("ISRES");
       LocalOptCombo->addItem("STOGO");

       LocalOptLayout->addWidget(LocalOptLabel);
       LocalOptLayout->addWidget(LocalOptCombo);

       // GNUPLOT
       QGridLayout *GNUplotLayout = new QGridLayout;
       QGroupBox *groupBoxGNUplot = new QGroupBox(tr("GNUplot settings"));
       GNUplotButton = new QPushButton("Browse");
       connect(GNUplotButton, SIGNAL(clicked()), this, SLOT(GNUplotOutput_clicked()));

       GNUplotLayout->addWidget(new QLabel("GNUplot output path:"), 0, 0);
       GNUplotLayout->addWidget(GNUplotButton, 0, 1);

       groupBoxGNUplot->setLayout(GNUplotLayout);


       QHBoxLayout * QucsLayout =  new QHBoxLayout();
       QucsLayout->addWidget(new QLabel("Qucs version"));
       QucsCombo = new QComboBox();
       QucsCombo->addItem("0.0.19");
       QucsCombo->addItem("0.0.18");
       QucsCombo->setEditable(true);
       QucsLayout->addWidget(QucsCombo);


      //Create go/cancel buttons
      RunButton = new QPushButton("Go");
      CancelButton = new QPushButton("Cancel");
      ButtonsLayout->addWidget(RunButton);
      ButtonsLayout->addWidget(CancelButton);

      connect(RunButton, SIGNAL(clicked()), this, SLOT(go_clicked()));
      connect(CancelButton, SIGNAL(clicked()), this, SLOT(cancel_clicked()));


      //Create main layout and add individual layouts
      QVBoxLayout * mainLayout = new QVBoxLayout();
      mainLayout->addLayout(Impedancelayout);
      mainLayout->addWidget(groupBox);
      mainLayout->addLayout(helpbox);
      mainLayout->addLayout(TopoLayout);
      mainLayout->addLayout(LocalOptLayout);
      mainLayout->addWidget(groupBoxGNUplot);
      mainLayout->addLayout(QucsLayout);
      mainLayout->addLayout(ButtonsLayout);

      centralWidget->setLayout(mainLayout);
      setCentralWidget(centralWidget);
      setWindowTitle("GRABIM");

      statusBar()->showMessage(tr("Ready"));
}

void ui::ArbitraryTopology_clicked()
{
    if (ArbitraryTopology->isChecked())
    {
        ArbitraryTopologyLineEdit->setEnabled(true);
        CodeLabel->setVisible(true);
        setFixedSize(400, 650);
    }
    else
    {
        ArbitraryTopologyLineEdit->setEnabled(false);
        CodeLabel->setVisible(false);
        setFixedSize(400, 500);
    }
}

void ui::SearchBestTopology_clicked()
{
        if (SearchBestTopology->isChecked())
        {
            ArbitraryTopologyLineEdit->setEnabled(false);
            CodeLabel->setVisible(false);
            setFixedSize(400, 450);
        }
        else
        {
            ArbitraryTopologyLineEdit->setEnabled(true);
            CodeLabel->setVisible(true);
            setFixedSize(400, 600);
        }

}

void ui::go_clicked()
{
    QMessageBox::information(0, QObject::tr("Started"),
                         QObject::tr("GRABIM has just started running. Please hang out with the terminal"));


    //Before starting the matching engine, we must ensure that the impedance data is already loaded
    if (SourceFile.isEmpty() && (FixedZSLineedit->text().isEmpty()))
    {
        QMessageBox::warning(0, QObject::tr("Error"),
                             QObject::tr("Please select a S-parameter file for the source termination"));
        return;
    }

    if (LoadFile.isEmpty()&& (FixedZLLineedit->text().isEmpty()))
    {
        QMessageBox::warning(0, QObject::tr("Error"),
                             QObject::tr("Please select a S-parameter file for the load termination"));
        return;
    }

    //Check whether the impedances are specified as s1p or s2p
    int formatSource=-1, formatLoad=-1;;
    if (SourceFile.contains(".s1p")) formatSource = 0;
    if (SourceFile.contains(".z1p")) formatSource = 1;
    if (SourceFile.contains(".s2p")) formatSource = 2;
    if (LoadFile.contains(".s1p")) formatLoad = 0;
    if (LoadFile.contains(".z1p")) formatLoad = 1;
    if (LoadFile.contains(".s2p")) formatLoad = 2;

    if ((formatSource != 0) && (formatSource != 1) && (FixedZSLineedit->text().isEmpty()))
    {
        QMessageBox::warning(0, QObject::tr("Error"),
                             QObject::tr("The source termination must be either a s1p or s2p file"));
        return;
    }

    if ((formatLoad != 0) && (formatLoad != 1) && (FixedZLLineedit->text().isEmpty()))
    {
        QMessageBox::warning(0, QObject::tr("Error"),
                             QObject::tr("The load termination must be either a s1p or s2p file"));
        return;
    }

    //Impedance data paths were already specified, let's proceed to bring the S-par data into memory
    IO inout_operations;

    if (!FixedZSCheckbox->isChecked())//Read source impedance from file
    {
       switch (formatSource)
       {
          case 0: inout_operations.loadS1Pdata(SourceFile.toStdString(), SOURCE, false); break;//s1p
          case 1: inout_operations.loadS1Pdata(SourceFile.toStdString(), SOURCE, true); break;//z1p
          case 2: inout_operations.loadS2Pdata(SourceFile.toStdString(), SOURCE);break;//s2p
       }
    }
    else//Set constant source impedance
    {
       cx_double zs_temp;
       QByteArray ba = FixedZSLineedit->text().toLatin1();
       char * text = ba.data();
       zs_temp = getComplexImpedanceFromText(text);
       if (zs_temp.real() == -1)//Check if the input value is correct
       {
           QMessageBox::warning(0, QObject::tr("Error"),
                                QObject::tr("The input given for the source impedance is not valid"));
           return;
       }
       inout_operations.set_constant_ZS_vs_freq(zs_temp);
    }

    if (!FixedZLCheckbox->isChecked())
    {
        switch (formatLoad)
        {
           case 0: inout_operations.loadS1Pdata(LoadFile.toStdString(), LOAD, false); break;//s1p
           case 1: inout_operations.loadS1Pdata(LoadFile.toStdString(), LOAD, true); break;//z1p
           case 2: inout_operations.loadS2Pdata(LoadFile.toStdString(), LOAD);break;//s2p
        }
    }
    else
    {
        cx_double zl_temp;
        QByteArray ba = FixedZLLineedit->text().toLatin1();
        char * text = ba.data();
        zl_temp = getComplexImpedanceFromText(text);

        if (zl_temp.real() == -1)//Check if the input value is correct
        {
            QMessageBox::warning(0, QObject::tr("Error"),
                                 QObject::tr("The input given for the load impedance is not valid"));
            return;
        }
        inout_operations.set_constant_ZL_vs_freq(zl_temp);

    }

    //Check frequency specifications
    double fmatching_min = minFEdit->text().toDouble();
    double fmatching_max = maxFEdit->text().toDouble();
    if ((fmatching_min == -1) || (fmatching_max == -1))
    {
        QMessageBox::warning(0, QObject::tr("Error"),
                             QObject::tr("Incorrect frequency settings"));
        return;
    }
    else//Everything correct... lets set frequency
    {
        //Scale frequency according to the combobox units
        fmatching_min *= getFreqScale(minFUnitsCombo->currentIndex());
        fmatching_max *= getFreqScale(maxFUnitsCombo->currentIndex());
        inout_operations.set_matching_band(fmatching_min, fmatching_max);

        //Check if the specified frequencies lie with the s1p/s2p data
        inout_operations.ResampleImpedances();//Force data update
        if (fmatching_min < inout_operations.getFrequency().min())//The lower freq is not present at s1p/s2p
        {
            QMessageBox::warning(0, QObject::tr("Error"),
                                 QObject::tr("One of the impedance data files does not contain the specified lower frequency"));
            return;
        }
        if (fmatching_max > inout_operations.getFrequency().max())//The maximum freq is not present at s1p/s2p
        {
            QMessageBox::warning(0, QObject::tr("Error"),
                                 QObject::tr("One of the impedance data files does not contain the specified upper frequency"));
            return;
        }
    }


    GRABIM MatchingObject;
    // Impedance and frequency settings
    MatchingObject.SetSourceImpedance(inout_operations.getSourceImpedance());
    MatchingObject.SetLoadImpedance(inout_operations.getLoadImpedance());
    MatchingObject.SetFrequency(inout_operations.getFrequency());

    //Topology
    if (ArbitraryTopology->isChecked())//Use custom topology
    {
        MatchingObject.SetTopology(ArbitraryTopologyLineEdit->text().toStdString());
    }
    else
    {//Search over the predefined circuit set
       MatchingObject.SetTopology("-1");
    }

    // Set local optimiser
    switch(LocalOptCombo->currentIndex())
    {
       case 0: inout_operations.setLocalOptimiser(nlopt::LN_NELDERMEAD);break;
       case 1: inout_operations.setLocalOptimiser(nlopt::LN_PRAXIS);break;
       case 2: inout_operations.setLocalOptimiser(nlopt::LN_SBPLX);break;
       case 3: inout_operations.setLocalOptimiser(nlopt::LN_COBYLA);break;
       case 4: inout_operations.setLocalOptimiser(nlopt::LN_BOBYQA);break;
       case 5: inout_operations.setLocalOptimiser(nlopt::LN_AUGLAG);break;
       case 6: inout_operations.setLocalOptimiser(nlopt::GN_ESCH);break;
       case 7: inout_operations.setLocalOptimiser(nlopt::GN_ISRES);break;
       case 8: inout_operations.setLocalOptimiser(nlopt::GD_STOGO);break;
    }



    GRABIM_Result R = MatchingObject.RunGRABIM();//Runs GRABIM. Well, this is not exactly the algorithm
    // detailed at [1] but from my point of view is functional and easier to code...
    //Notes:
    // 1) The candidate vector is not in log units. I do not see a good reason for doing so. Maybe I am missing something important
    // 2) Frequency is not in rad/s.
    // 3) The objective function is the magnitude of S11 expressed in dB. log(x) functions usually have strong
    // gradients so it seem to suggest that this is good for derivative free opt algorithms
    // 4) This code takes advantage from NLopt derivative-free local optimisers. NLopt is easy to link and it works
    // fine. Despite the fact that the Nelder-Mead algorithm does not guarantee convergence (among other problems), it leads to achieve a good local
    // (probably, global) optimum. This is caused by the fact that the matching network should be as simple as possible => few elements => xk \in R^N, where
    // N is typically < 6. Even N=6 is a big number, please consider that matching networks are tight to physical constraints in practice, so, the larger the
    // network, the harder the 'tuning'.
    R.source_path = SourceFile;
    R.load_path = LoadFile;
    R.QucsVersion = QucsCombo->currentText();



    if (GNUplot_path.isEmpty())
    {
        GNUplot_path = QCoreApplication::applicationDirPath() + "GRABIM.dat";
    }

    QMessageBox::information(0, QObject::tr("Finished"),
                         QObject::tr("GRABIM has successfully finished. \nPlease execute: 'gnuplot plotscript' to take a look to the results.\nThe network has been copied into the clipboard as a Qucs netlist"));

    inout_operations.exportGNUplot(R, GNUplot_path.toStdString());
    inout_operations.ExportQucsSchematic(R);
}

void ui::cancel_clicked()
{
   this->close();
}

void ui::SourceImpedance_clicked()
{
    SourceFile = QFileDialog::getOpenFileName(this,
        tr("Open S-parameter file"), QCoreApplication::applicationDirPath());
}
void ui::LoadImpedance_clicked()
{
    LoadFile = QFileDialog::getOpenFileName(this,
        tr("Open S-parameter file"), QCoreApplication::applicationDirPath());
}

// GNUplot output path
void ui::GNUplotOutput_clicked()
{
    GNUplot_path = QFileDialog::getSaveFileName(this,
        tr("Open S-parameter file"), QCoreApplication::applicationDirPath());
}

// Shows hides button/lineedit for source impedance selection
void ui::FixedZSCheckbox_clicked()
{
    if (FixedZSCheckbox->isChecked())
    {
        FixedZSLineedit->setVisible(true);
        ZSOhmLabel->setVisible(true);
        SourceFileButton->setVisible(false);
    }
    else
    {
        FixedZSLineedit->setVisible(false);
        ZSOhmLabel->setVisible(false);
        SourceFileButton->setVisible(true);
    }
}

// Shows hides button/lineedit for load impedance selection
void ui::FixedZLCheckbox_clicked()
{
    if (FixedZLCheckbox->isChecked())
    {
        FixedZLLineedit->setVisible(true);
        ZLOhmLabel->setVisible(true);
        LoadFileButton->setVisible(false);
    }
    else
    {
        FixedZLLineedit->setVisible(false);
        ZLOhmLabel->setVisible(false);
        LoadFileButton->setVisible(true);
    }
}

// The user can also specify a constant complex impedance vs frequency. This function parses user's input
// and return the cx_double value
cx_double ui::getComplexImpedanceFromText(char *arg)
{
    string arg_str = arg;
    int index = arg_str.find_first_of("j");
    int sign = 1;
    double zreal, zimag;
    if (index != -1)//Then it is a single impedance
    {
        zreal = atof(arg_str.substr(0, index-1).c_str());
        if (!arg_str.substr(index-1, 1).compare("-")) sign = -1;
        zimag = sign*atof(arg_str.substr(index+1).c_str());
        return cx_double(zreal, zimag);
    }
    else
    {
        zreal = atof(arg);
        if (zreal > 0)return cx_double(zreal, 0);
        else
            return cx_double(-1, -1);
    }
}

// Gets frequency scale from the combobox index
double ui::getFreqScale(int index)
{
   switch (index)
   {
      case 0: return 1e3;break;
      case 1: return 1e6;break;
      case 2: return 1e9;break;

   }
   return 1e6;
}



