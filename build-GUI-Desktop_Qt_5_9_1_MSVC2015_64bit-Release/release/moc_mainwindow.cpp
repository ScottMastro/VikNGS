/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../RVS/mainwindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#include <QtCore/QVector>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_Runner_t {
    QByteArrayData data[13];
    char stringdata0[111];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Runner_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Runner_t qt_meta_stringdata_Runner = {
    {
QT_MOC_LITERAL(0, 0, 6), // "Runner"
QT_MOC_LITERAL(1, 7, 8), // "complete"
QT_MOC_LITERAL(2, 16, 0), // ""
QT_MOC_LITERAL(3, 17, 7), // "process"
QT_MOC_LITERAL(4, 25, 6), // "replot"
QT_MOC_LITERAL(5, 32, 19), // "std::vector<double>"
QT_MOC_LITERAL(6, 52, 6), // "values"
QT_MOC_LITERAL(7, 59, 5), // "setUi"
QT_MOC_LITERAL(8, 65, 15), // "Ui::MainWindow*"
QT_MOC_LITERAL(9, 81, 2), // "ui"
QT_MOC_LITERAL(10, 84, 10), // "setRequest"
QT_MOC_LITERAL(11, 95, 7), // "Request"
QT_MOC_LITERAL(12, 103, 7) // "request"

    },
    "Runner\0complete\0\0process\0replot\0"
    "std::vector<double>\0values\0setUi\0"
    "Ui::MainWindow*\0ui\0setRequest\0Request\0"
    "request"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Runner[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   39,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       3,    0,   40,    2, 0x0a /* Public */,
       4,    1,   41,    2, 0x0a /* Public */,
       7,    1,   44,    2, 0x0a /* Public */,
      10,    1,   47,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 5,    6,
    QMetaType::Void, 0x80000000 | 8,    9,
    QMetaType::Void, 0x80000000 | 11,   12,

       0        // eod
};

void Runner::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Runner *_t = static_cast<Runner *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->complete(); break;
        case 1: _t->process(); break;
        case 2: _t->replot((*reinterpret_cast< std::vector<double>(*)>(_a[1]))); break;
        case 3: _t->setUi((*reinterpret_cast< Ui::MainWindow*(*)>(_a[1]))); break;
        case 4: _t->setRequest((*reinterpret_cast< Request(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (Runner::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&Runner::complete)) {
                *result = 0;
                return;
            }
        }
    }
}

const QMetaObject Runner::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_Runner.data,
      qt_meta_data_Runner,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *Runner::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Runner::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_Runner.stringdata0))
        return static_cast<void*>(const_cast< Runner*>(this));
    return QObject::qt_metacast(_clname);
}

int Runner::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 5)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 5;
    }
    return _id;
}

// SIGNAL 0
void Runner::complete()
{
    QMetaObject::activate(this, &staticMetaObject, 0, nullptr);
}
struct qt_meta_stringdata_MainWindow_t {
    QByteArrayData data[32];
    char stringdata0[622];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MainWindow_t qt_meta_stringdata_MainWindow = {
    {
QT_MOC_LITERAL(0, 0, 10), // "MainWindow"
QT_MOC_LITERAL(1, 11, 14), // "calculatePower"
QT_MOC_LITERAL(2, 26, 0), // ""
QT_MOC_LITERAL(3, 27, 15), // "QVector<double>"
QT_MOC_LITERAL(4, 43, 4), // "pval"
QT_MOC_LITERAL(5, 48, 10), // "sim_replot"
QT_MOC_LITERAL(6, 59, 6), // "values"
QT_MOC_LITERAL(7, 66, 4), // "maxX"
QT_MOC_LITERAL(8, 71, 15), // "constructGroups"
QT_MOC_LITERAL(9, 87, 35), // "std::vector<SimulationRequest..."
QT_MOC_LITERAL(10, 123, 3), // "run"
QT_MOC_LITERAL(11, 127, 6), // "ntests"
QT_MOC_LITERAL(12, 134, 16), // "constructRequest"
QT_MOC_LITERAL(13, 151, 17), // "SimulationRequest"
QT_MOC_LITERAL(14, 169, 6), // "groups"
QT_MOC_LITERAL(15, 176, 21), // "on_sim_runBtn_clicked"
QT_MOC_LITERAL(16, 198, 26), // "on_sim_groupAddBtn_clicked"
QT_MOC_LITERAL(17, 225, 29), // "on_sim_groupRemoveBtn_clicked"
QT_MOC_LITERAL(18, 255, 31), // "on_sim_testBootChk_stateChanged"
QT_MOC_LITERAL(19, 287, 4), // "arg1"
QT_MOC_LITERAL(20, 292, 26), // "on_sim_testBootChk_toggled"
QT_MOC_LITERAL(21, 319, 7), // "checked"
QT_MOC_LITERAL(22, 327, 30), // "on_sim_testRareCastBtn_toggled"
QT_MOC_LITERAL(23, 358, 32), // "on_sim_testRareCalphaBtn_toggled"
QT_MOC_LITERAL(24, 391, 25), // "on_main_vcfDirBtn_clicked"
QT_MOC_LITERAL(25, 417, 28), // "on_main_sampleDirBtn_clicked"
QT_MOC_LITERAL(26, 446, 25), // "on_main_bedDirBtn_clicked"
QT_MOC_LITERAL(27, 472, 22), // "on_main_runBtn_clicked"
QT_MOC_LITERAL(28, 495, 31), // "on_main_testRareCastBtn_toggled"
QT_MOC_LITERAL(29, 527, 33), // "on_main_testRareCalphaBtn_tog..."
QT_MOC_LITERAL(30, 561, 32), // "on_main_testBootChk_stateChanged"
QT_MOC_LITERAL(31, 594, 27) // "on_main_testBootChk_toggled"

    },
    "MainWindow\0calculatePower\0\0QVector<double>\0"
    "pval\0sim_replot\0values\0maxX\0constructGroups\0"
    "std::vector<SimulationRequestGroup>\0"
    "run\0ntests\0constructRequest\0"
    "SimulationRequest\0groups\0on_sim_runBtn_clicked\0"
    "on_sim_groupAddBtn_clicked\0"
    "on_sim_groupRemoveBtn_clicked\0"
    "on_sim_testBootChk_stateChanged\0arg1\0"
    "on_sim_testBootChk_toggled\0checked\0"
    "on_sim_testRareCastBtn_toggled\0"
    "on_sim_testRareCalphaBtn_toggled\0"
    "on_main_vcfDirBtn_clicked\0"
    "on_main_sampleDirBtn_clicked\0"
    "on_main_bedDirBtn_clicked\0"
    "on_main_runBtn_clicked\0"
    "on_main_testRareCastBtn_toggled\0"
    "on_main_testRareCalphaBtn_toggled\0"
    "on_main_testBootChk_stateChanged\0"
    "on_main_testBootChk_toggled"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MainWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      19,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,  109,    2, 0x08 /* Private */,
       5,    2,  112,    2, 0x08 /* Private */,
       8,    2,  117,    2, 0x08 /* Private */,
      12,    1,  122,    2, 0x08 /* Private */,
      15,    0,  125,    2, 0x08 /* Private */,
      16,    0,  126,    2, 0x08 /* Private */,
      17,    0,  127,    2, 0x08 /* Private */,
      18,    1,  128,    2, 0x08 /* Private */,
      20,    1,  131,    2, 0x08 /* Private */,
      22,    1,  134,    2, 0x08 /* Private */,
      23,    1,  137,    2, 0x08 /* Private */,
      24,    0,  140,    2, 0x08 /* Private */,
      25,    0,  141,    2, 0x08 /* Private */,
      26,    0,  142,    2, 0x08 /* Private */,
      27,    0,  143,    2, 0x08 /* Private */,
      28,    1,  144,    2, 0x08 /* Private */,
      29,    1,  147,    2, 0x08 /* Private */,
      30,    1,  150,    2, 0x08 /* Private */,
      31,    1,  153,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Double, 0x80000000 | 3,    4,
    QMetaType::Void, 0x80000000 | 3, QMetaType::Int,    6,    7,
    0x80000000 | 9, QMetaType::Int, QMetaType::Int,   10,   11,
    0x80000000 | 13, 0x80000000 | 9,   14,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   19,
    QMetaType::Void, QMetaType::Bool,   21,
    QMetaType::Void, QMetaType::Bool,   21,
    QMetaType::Void, QMetaType::Bool,   21,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,   21,
    QMetaType::Void, QMetaType::Bool,   21,
    QMetaType::Void, QMetaType::Int,   19,
    QMetaType::Void, QMetaType::Bool,   21,

       0        // eod
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        MainWindow *_t = static_cast<MainWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: { double _r = _t->calculatePower((*reinterpret_cast< QVector<double>(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< double*>(_a[0]) = std::move(_r); }  break;
        case 1: _t->sim_replot((*reinterpret_cast< QVector<double>(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 2: { std::vector<SimulationRequestGroup> _r = _t->constructGroups((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])));
            if (_a[0]) *reinterpret_cast< std::vector<SimulationRequestGroup>*>(_a[0]) = std::move(_r); }  break;
        case 3: { SimulationRequest _r = _t->constructRequest((*reinterpret_cast< std::vector<SimulationRequestGroup>(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< SimulationRequest*>(_a[0]) = std::move(_r); }  break;
        case 4: _t->on_sim_runBtn_clicked(); break;
        case 5: _t->on_sim_groupAddBtn_clicked(); break;
        case 6: _t->on_sim_groupRemoveBtn_clicked(); break;
        case 7: _t->on_sim_testBootChk_stateChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 8: _t->on_sim_testBootChk_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 9: _t->on_sim_testRareCastBtn_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 10: _t->on_sim_testRareCalphaBtn_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 11: _t->on_main_vcfDirBtn_clicked(); break;
        case 12: _t->on_main_sampleDirBtn_clicked(); break;
        case 13: _t->on_main_bedDirBtn_clicked(); break;
        case 14: _t->on_main_runBtn_clicked(); break;
        case 15: _t->on_main_testRareCastBtn_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 16: _t->on_main_testRareCalphaBtn_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 17: _t->on_main_testBootChk_stateChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 18: _t->on_main_testBootChk_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 0:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QVector<double> >(); break;
            }
            break;
        case 1:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QVector<double> >(); break;
            }
            break;
        }
    }
}

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow.data,
      qt_meta_data_MainWindow,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow.stringdata0))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 19)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 19;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 19)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 19;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
