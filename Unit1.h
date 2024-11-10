
#ifndef Unit1H
#define Unit1H

#define VC_EXTRALEAN
#define WIN32_EXTRA_LEAN
#define WIN32_LEAN_AND_MEAN

#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <ComCtrls.hpp>
#include <Buttons.hpp>
#include <Menus.hpp>
#include <Dialogs.hpp>
#include "CSPIN.h"

#include <vector>
#include <stdint.h>

#define WM_INIT_GUI						(WM_USER + 100)
#define WM_COMPUTE_DISH 				(WM_USER + 101)

#define ARRAY_SIZE(array)           (sizeof(array) / sizeof(array[0]))

typedef struct t_dish_point
{
	t_dish_point()
   {
   	r = 0.0;
      y = 0.0;
   }
	t_dish_point(const double _r, const double _y)
	{
   	r = _r;
      y = _y;
   }
	double r;   // radius from center
	double y;   // rear dish parabaloid curve Y
} t_dish_point;

typedef struct t_petal_point
{
	t_petal_point()
   {
   	x = 0.0;
      y = 0.0;
   	r = 0.0;
      s = 0.0;
      w = 0.0;
   }
	t_petal_point(const double _x, const double _y, const double _r, const double _s, const double _w)
	{
   	x = _x;
      y = _y;
   	r = _r;
      s = _s;
      w = _w;
   }
	double x;   // rear dish parabaloid curve radius from center
	double y;   // rear dish parabaloid curve Y
	double r;   // radius from center
   double s;   // petal cut out width
	double w;   // petal segment width
} t_petal_point;

class TForm1 : public TForm
{
__published:	// IDE-managed Components
	TOpenDialog *OpenDialog1;
	TSaveDialog *SaveDialog1;
	TPaintBox *PaintBox1;
	TPanel *Panel1;
	TMemo *Memo1;
	TPanel *Panel2;
	TEdit *FrequencyEdit;
	TLabel *Label3;
	TLabel *Label4;
	TLabel *Label6;
	TLabel *Label7;
	TLabel *Label8;
	TLabel *Label9;
	TEdit *FDRatioEdit;
	TLabel *Label11;
	TLabel *Label1;
	TEdit *DiameterEdit;
	TLabel *Label2;
	TCSpinEdit *PointsCSpinEdit;
	TCSpinEdit *PetalsCSpinEdit;
	TCSpinEdit *EfficiencyCSpinEdit;
	void __fastcall FormCreate(TObject *Sender);
	void __fastcall FormDestroy(TObject *Sender);
	void __fastcall FormClose(TObject *Sender, TCloseAction &Action);
	void __fastcall PaintBox1Paint(TObject *Sender);
	void __fastcall FormKeyDown(TObject *Sender, WORD &Key,
			 TShiftState Shift);
	void __fastcall PaintBox1MouseMove(TObject *Sender, TShiftState Shift,
          int X, int Y);
	void __fastcall PaintBox1MouseDown(TObject *Sender, TMouseButton Button,
          TShiftState Shift, int X, int Y);
	void __fastcall PaintBox1MouseUp(TObject *Sender, TMouseButton Button,
			 TShiftState Shift, int X, int Y);
	void __fastcall EditKeyDown(TObject *Sender, WORD &Key,
          TShiftState Shift);
	void __fastcall EditKeyUp(TObject *Sender, WORD &Key,
          TShiftState Shift);
	void __fastcall CSpinChange(TObject *Sender);

private:
	String             m_ini_filename;

	int                m_screen_width;
	int                m_screen_height;
	SYSTEM_INFO        m_system_info;

	Graphics::TBitmap *m_graph_bm;

	double             m_frequency_Hz;
	double             m_gain_dBi;
	double             m_diameter_meters;
	double             m_focus_diameter_ratio;
	double             m_efficiency;
	int                m_petals;
	int                m_points;

	double             m_area_sq_meters;
	double             m_gain;
	double             m_wave_length_meters;
	double             m_focus_point_meters;
	double             m_gain_absolute;
	double             m_feed_angle;
	double             m_feed_arm_meters;
	double             m_space_attenuation_dB;
	double             m_desired_taper_dB;
	double             m_beam_width_3dB;
	double             m_simple_horn_diameter;

	double             m_waveguide_h;
	double             m_waveguide_e;
	double             m_horn_gain_dBi;

	std::vector <t_dish_point>  m_dish_point;
	std::vector <t_petal_point> m_petal_point;	// dimesions for a flat sheet cut out
	std::vector < std::pair <double, double> > m_beam_width;

	int                m_mouse_x;
	int                m_mouse_y;

	float              m_pz_offset_x;
	float              m_pz_offset_y;
	float              m_pz_scale_x;
	float              m_pz_scale_y;

	int                m_left_margin;
	int                m_right_margin;

	void __fastcall loadSettings();
	void __fastcall saveSettings();

	void __fastcall comboBoxAutoWidth(TComboBox *comboBox);

	void __fastcall doUpdate();

	void __fastcall drawArc(Graphics::TBitmap *bm, const int x, const int y, const double start_deg, const double stop_deg, const double radius);
	void __fastcall drawCircle(Graphics::TBitmap *bm, const int x, const int y, const int radius);

	void __fastcall WMWindowPosChanging(TWMWindowPosChanging &msg);
	void __fastcall CMMouseEnter(TMessage &msg);
	void __fastcall CMMouseLeave(TMessage &msg);
	void __fastcall WMInitGUI(TMessage &msg);
	void __fastcall WMComputeDish(TMessage &msg);

protected:
	#pragma option push -vi-
	BEGIN_MESSAGE_MAP
		VCL_MESSAGE_HANDLER(WM_WINDOWPOSCHANGING, TWMWindowPosMsg, WMWindowPosChanging);

		VCL_MESSAGE_HANDLER(CM_MOUSELEAVE, TMessage, CMMouseLeave);
		VCL_MESSAGE_HANDLER(CM_MOUSEENTER, TMessage, CMMouseEnter);

		VCL_MESSAGE_HANDLER(WM_INIT_GUI, TMessage, WMInitGUI);

		VCL_MESSAGE_HANDLER(WM_COMPUTE_DISH, TMessage, WMComputeDish);

	END_MESSAGE_MAP(TForm)
	#pragma option pop

public:		// User declarations
	__fastcall TForm1(TComponent* Owner);
};

extern PACKAGE TForm1 *Form1;

#endif
