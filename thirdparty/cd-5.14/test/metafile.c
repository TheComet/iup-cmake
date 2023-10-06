#include <cd.h>
#include <cdmf.h>
#include <cdcgm.h>

void draw1(cdCanvas *canvas)
{
	cdCanvasMarkSize(canvas, 5);
	cdCanvasMarkType(canvas, CD_PLUS);
	cdCanvasMark(canvas, 10,90);
	cdCanvasMarkType(canvas, CD_STAR);
	cdCanvasMark(canvas, 20,90);
	cdCanvasMarkType(canvas, CD_CIRCLE);
	cdCanvasMark(canvas, 30,90);
	cdCanvasMarkType(canvas, CD_X);
	cdCanvasMark(canvas, 40,90);
	cdCanvasMarkType(canvas, CD_BOX);
	cdCanvasMark(canvas, 50,90);
	cdCanvasMarkType(canvas, CD_DIAMOND);
	cdCanvasMark(canvas, 60,90);
	cdCanvasMarkType(canvas, CD_HOLLOW_CIRCLE);
	cdCanvasMark(canvas, 70,90);
	cdCanvasMarkType(canvas, CD_HOLLOW_BOX);
	cdCanvasMark(canvas, 80,90);
	cdCanvasMarkType(canvas, CD_HOLLOW_DIAMOND);
	cdCanvasMark(canvas, 90,90);

	cdCanvasLineStyle(canvas, CD_CONTINUOUS);
	cdCanvasLine(canvas, 10,80,80,80);
	cdCanvasLineStyle(canvas, CD_DASHED);
	cdCanvasLine(canvas, 10,75,80,75);
	cdCanvasLineStyle(canvas, CD_DOTTED);
	cdCanvasLine(canvas, 10,70,80,70);
	cdCanvasLineStyle(canvas, CD_DASH_DOT);
	cdCanvasLine(canvas, 10,65,80,65);
	cdCanvasLineStyle(canvas, CD_DASH_DOT_DOT);
	cdCanvasLine(canvas, 10,60,80,60);

	cdCanvasLineStyle(canvas, CD_CONTINUOUS);

	cdCanvasHatch(canvas, CD_HORIZONTAL);

	cdCanvasBegin(canvas, CD_FILL);
	cdCanvasVertex(canvas, 10,50);
	cdCanvasVertex(canvas, 50,50);
	cdCanvasVertex(canvas, 50,10);
	cdCanvasVertex(canvas, 10,10);
	cdCanvasEnd(canvas);

	cdCanvasHatch(canvas, CD_VERTICAL);
	cdCanvasBegin(canvas, CD_FILL);
	cdCanvasVertex(canvas, 60,50);
	cdCanvasVertex(canvas, 100,50);
	cdCanvasVertex(canvas, 100,10);
	cdCanvasVertex(canvas, 60,10);
	cdCanvasEnd(canvas);	

	cdCanvasHatch(canvas, CD_FDIAGONAL);
	cdCanvasBegin(canvas, CD_FILL);
	cdCanvasVertex(canvas, 110,50);
	cdCanvasVertex(canvas, 150,50);
	cdCanvasVertex(canvas, 150,10);
	cdCanvasVertex(canvas, 110,10);
	cdCanvasEnd(canvas);	

	cdCanvasHatch(canvas, CD_BDIAGONAL);
	cdCanvasBegin(canvas, CD_FILL);
	cdCanvasVertex(canvas, 160,50);
	cdCanvasVertex(canvas, 200,50);
	cdCanvasVertex(canvas, 200,10);
	cdCanvasVertex(canvas, 160,10);
	cdCanvasEnd(canvas);	

	cdCanvasHatch(canvas, CD_CROSS);
	cdCanvasBegin(canvas, CD_FILL);
	cdCanvasVertex(canvas, 210,50);
	cdCanvasVertex(canvas, 250,50);
	cdCanvasVertex(canvas, 250,10);
	cdCanvasVertex(canvas, 210,10);
	cdCanvasEnd(canvas);	

	cdCanvasHatch(canvas, CD_DIAGCROSS);
	cdCanvasBegin(canvas, CD_FILL);
	cdCanvasVertex(canvas, 260,50);
	cdCanvasVertex(canvas, 300,50);
	cdCanvasVertex(canvas, 300,10);
	cdCanvasVertex(canvas, 260,10);
	cdCanvasEnd(canvas);	
	
	cdCanvasFont(canvas, "System",CD_BOLD,CD_STANDARD);
	cdCanvasText(canvas, 10,100,"Text");
	cdCanvasFont(canvas, "Courier",CD_BOLD,CD_STANDARD);
	cdCanvasText(canvas, 60,100,"Text");
	cdCanvasFont(canvas, "Times",CD_BOLD,CD_STANDARD);
	cdCanvasText(canvas, 110,100,"Text");
	cdCanvasFont(canvas, "Helvetica",CD_BOLD,CD_STANDARD);
	cdCanvasText(canvas, 160,100,"Text");
}

void main(void)
{
//	cdCanvas *canvas = cdCreateCanvas(CD_METAFILE,"TESTE.MF 100x100");
	cdCanvas *canvas = cdCreateCanvas(CD_CGM,"simple1.cgm 100x100 5 -t");  /* 500x500 pixels */
	draw1(canvas);
	cdKillCanvas(canvas);
}


