import ROOT
from PyQt6.QtCore import QDateTime, Qt, QTimer
from PyQt6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDateTimeEdit,
        QDial, QDialog, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
        QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QSpinBox, QStyleFactory, QTableWidget, QTabWidget, QTextEdit,
        QVBoxLayout, QWidget,QMainWindow)

import sys

def tree_columns():
    
    PATH = "/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/49Ti_analyzed/"
    treeName = "SPSTree"
    
    Globbed_Raw_Frame = ROOT.RDataFrame(treeName, PATH + "*.root")
    raw_columns = Globbed_Raw_Frame.GetColumnNames()
    print(raw_columns)
    
            
app = QApplication(sys.argv)

raw_columns_button = QPushButton("Print Columns in SPSTree")
raw_columns_button.clicked.connect(tree_columns)
raw_columns_button.show()

app.exec()
    