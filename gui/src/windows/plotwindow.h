#ifndef PLOTWINDOW_H
#define PLOTWINDOW_H

#include <QWidget>
#include "../src/Variant.h"
#include "../widgets/qcustomplot.h"

class Chromosome{
private:
    QVector<Variant> variants;
    QColor colour;
    QString name;
    double maxPos;
    double minPos;
    double offset;

public:
    Chromosome(){
        name="unknown";
        maxPos = -1;
        minPos = 2100000000;
    }

    Chromosome(QString name, Variant first, QColor c){
        this->name = name;
        variants.push_back(first);
        maxPos = first.pos;
        minPos = maxPos;
        colour = c;
    }

    void addVariant(Variant toAdd){
        variants.push_back(toAdd);
        if(maxPos < toAdd.pos)
            maxPos = toAdd.pos;
        if(minPos > toAdd.pos)
            minPos = toAdd.pos;
    }

    QVector<double> getPositions(double offset){

        QVector<double> pos;
        for(int i = 0; i < variants.size(); i++)
            pos.push_back(offset + variants[i].pos);

        return pos;
    }
    QVector<double> getPositions(){ return getPositions(0); }


    QVector<double> getPvals(){

        QVector<double> pvals;
        for(int i = 0; i < variants.size(); i++)
            pvals.push_back(-log10(variants[i].pvalue));

        return pvals;
    }

    Variant getVariant(int i){return variants[i];}
    void setColour(QColor c){ this->colour = c; }
    QColor getColour(){ return this->colour; }
    int size(){ return variants.size(); }
    double getSpan(){ return maxPos - minPos; }
    double getMaxPos(){ return maxPos; }
    QString getName(){ return name; }
    void setName(QString n){ this->name = n; }

};


namespace Ui {
class PlotWindow;
}

class PlotWindow : public QWidget
{
    Q_OBJECT

public:
    explicit PlotWindow(QWidget *parent = 0);
    ~PlotWindow();

public slots:

    void initialize(QVector<Variant> values, QString title);
    void createChromosomes(QVector<Variant> variants);

    //for development
    void setRandomChromosomes();
    Chromosome generateRandomChromosome(int n, std::string chrom, int maxPos);
    //

    void mouseMoveWindow(QMouseEvent* event);
    void mouseMoveGenome(QMouseEvent* event);
    void mouseClickGenome(QMouseEvent *event);
    QString getChromUnderCursor(QMouseEvent *event);
    void mouseMoveChromosome(QMouseEvent *event);
    void mouseScrollChromosome(QWheelEvent* event);
    void mouseClickChromosome(QMouseEvent *event);

    void buildGenomePlot();
    void zoomChromosomePlot(double min, double max);
    void buildChromosomePlot(QString chrName, double min, double max);
    void updateVariantHighlightLayer(Variant variants);

    void updateVariantInfo(Variant variant);

    QCPGraph* getGraphByName(QCustomPlot *plot, QString name);
    Variant findClosestVariant(double x, double y, double maxDist);

private:
    Ui::PlotWindow *ui;  
    QColor grey1 = QColor::fromRgb(190, 190, 190);
    QColor grey2 = QColor::fromRgb(169, 169, 169);
    QColor highlight = QColor::fromRgb(255, 127, 80);
    QColor focus = QColor::fromRgb(102, 204, 204);

    QMap<QString, Chromosome> chromosomes;
    QVector<QString> chrNames;
    QCPItemLine *horizontal;
    QString focusedChr;
    Variant focusedVar;
    Variant nullVariant;

    Variant highlightVar;

    const QString chrGraphName = "chromgraph";

    static int extractNumber(const QString &c)
    {
        QString number = "";

        for(int i=0; i < c.size(); i++)
            if(c[i].isDigit())
                number += c[i];

        if(number == "")
            return -1;
        return number.toInt();

    }

    static bool chrCompare(const QString &c1, const QString &c2)
    {
        int i1 = extractNumber(c1);
        int i2 = extractNumber(c2);

        if(i1 > 0 && i2 > 0)
            return i1 < i2;

        return c1.toLower() < c2.toLower();
    }

    bool eventFilter(QObject *obj, QEvent *event)
    {
      if (event->type() == QEvent::MouseMove)
      {
        QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
        mouseMoveWindow(mouseEvent);
      }
      return false;
    }
};

#endif // PLOTWINDOW_H