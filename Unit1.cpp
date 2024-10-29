
#include <vcl.h>
#include <inifiles.hpp>
#include <jpeg.hpp>

#include <stdio.h>
#include <fastmath.h>
#include <limits.h>
#include <stdlib.h>     // srand, rand

#pragma hdrstop

#include "Unit1.h"

#pragma package(smart_init)
#pragma link "CSPIN"
#pragma resource "*.dfm"

#if defined(__BORLANDC__) && (__BORLANDC__ < 0x0600)
	#define floorf	      (float)floor
	#define fmodf	      (float)fmod
	#define sqrtf        (float)sqrt
	#define sinf         (float)sin
	#define cosf         (float)cos
 	#define tanf         (float)tan
	#define atan2f       (float)atan2
	#define ceilf        (float)ceil
	#define log10f       (float)log10
	#define logf         (float)log
	#define powf         (float)pow
	#define expf         (float)exp
	#define fabsf        (float)fabs
	#define asinf        (float)asin
	#define acosf        (float)acos
#endif

#define max(a, b)       ((a > b) ? a : b)
#define min(a, b)       ((a < b) ? a : b)
//#define max(a, b)       ({__typeof__ (a) _a = (a); __typeof__ (b) _b = (b); (_a > _b) ? _a : _b; })
//#define min(a, b)       ({__typeof__ (a) _a = (a); __typeof__ (b) _b = (b); (_a < _b) ? _a : _b; })

const double sol = 299792458;   // meters per second

const double deg_to_rad = M_PI / 180.0;
const double rad_to_deg = 180.0 / M_PI;

const int MaxPoints = 1000;
const int MaxPetals = 1000;

TForm1 *Form1 = NULL;

const ERR_MARGIN = 2e-3;

typedef struct
{
	uint16_t MajorVer;
	uint16_t MinorVer;
	uint16_t ReleaseVer;
	uint16_t BuildVer;
} TVersion;

// ************************************************************************

std::vector <String> __fastcall stringSplit(String s, String separator)
{
	std::vector <String> strings;

	if (separator.IsEmpty())
		return strings;

	while (s.Length() > 0)
	{
		const int p = s.Pos(separator);
		if (p <= 0)
			break;
		String s2 = s.SubString(1, p - 1);
		s = s.SubString(p + separator.Length(), s.Length());
		strings.push_back(s2);
	}

	if (!s.IsEmpty())
		strings.push_back(s);

	return strings;
}

// dB to linear
__inline float __fastcall dB2lin(const float dB)
{
	return powf(10.0f, dB / 10.0f);
}

 // linear to dB
__inline float __fastcall lin2dB(const float lin)
{
	return 10.0f * log10f(lin);
}

double __fastcall phase_error(const double aperture, const double flare_angle, const double wavelength)
{
	const double rho = 0.5 * aperture / sin(flare_angle);
	return aperture * aperture / (8 * wavelength * rho);
}

// MathCad fit by Matt Reilly of Fig. 23 from Balanis
//
// input is Maximum input phase error in wavelengths
// returns loss due to phase error in dB.
//
//	Le(s) = sum(from i = 1 to 10) Xle[i] * s^i
double __fastcall Le(const double s)
{
	const double Xle[] = {
		10.5075396618,
		-236.610713851,
		2603.14290758,
		-15318.8133545,
		55318.7209697,
		-125471.932662,
		178426.752358,
		-153811.177005,
		73233.5757222,
		-14743.1657614
	};

	if (s < 0.0 || s > 1.0)
	{
		//cerr << "ERROR: excessive phase error in Le\n\n";
		return 1;
	}

	double le = 0.0;
	for (int i = 0; i < ARRAY_SIZE(Xle); i++)
		le += Xle[i] * pow(s, 1 + i);
		//le += Xle[i] * Power(s, 1 + i);

	return le;
}

// MathCad fit by Matt Reilly of Fig. 23 from Balanis
//
// input is Maximum input phase error in wavelengths
// returns loss due to phase error in dB.
//
//	Lh(t) = sum(from i=1 to 10) Xlh[i] * t^i
double __fastcall Lh(const double t)
{
	const double Xlh[] = {
		-13.9134920465,
		398.4347217,
		-4241.37951372,
		24010.7335428,
		-79636.863311,
		162957.754379,
		-208325.065791,
		162037.036743,
		-70133.3772858,
		12951.940007
	};

	if (t < 0.0 || t > 1.0)
	{
		//cerr << "ERROR: excessive phase error in Lh\n\n";
		return 1;
	}

	double lh = 0.0;
	for (int i = 0; i < ARRAY_SIZE(Xlh); i++)
		lh += Xlh[i] * pow(t, 1 + i);
		//lh += Xlh[i] * Power(t, 1 + i);

	return lh;
}

double __fastcall epc(const double l, const double q, const double c2, const double s2, const double c, const double s)
{
	return l * (1.0 - (q * ((c2 * c) + (s2 * s)) / ((c * c) + (s * s))));
}

double __fastcall hpc(const double l, const double r, const double t, const double c, const double s)
{
	return l * (1.0 + (((r * c) + (t * s)) / ((c * c) + (s * s))));
}

// from stegun and abramowitz pp 301 and 302.
double __fastcall fresnel_f(const double z)
{
	return (1.0 + (0.926 * z)) / (2.0 + (1.792 * z) + (3.104 * z * z));
}

double __fastcall fresnel_g(const double z)
{
	return 1.0 / (2.0 + (4.142 * z) + (3.492 * z * z) + (6.670 * z * z * z));
}

double __fastcall fresnel_cos(double z, double &low, double &high)
{
	double sign = 1;
	if (z < 0.0)
	{
		z = fabs(z);
		sign = -1;
	}

	const double x = (M_PI / 2) * z * z;
	const double f = fresnel_f(z);
	const double g = fresnel_g(z);
	const double c = cos(x);
	const double s = sin(x);

	const double sum = 0.5 + (f * s) - (g * c);

	const double fs1 = (f + ERR_MARGIN) * s;
	const double fs2 = (f - ERR_MARGIN) * s;
	const double gc1 = (g + ERR_MARGIN) * c;
	const double gc2 = (g - ERR_MARGIN) * c;

	if (sign > 0)
	{
		low  =  (0.5 + min(fs1, fs2) - max(gc1, gc2));
		high =  (0.5 + max(fs1, fs2) - min(gc1, gc2));
	}
	else
	{
		low  = -(0.5 + max(fs1, fs2) - min(gc1, gc2));
		high = -(0.5 + min(fs1, fs2) - max(gc1, gc2));
	}

	return sum * sign;
}

double __fastcall fresnel_sin(double z, double &low, double &high)
{
	double sign = 1;
	if (z < 0.0)
	{
		z = fabs(z);
		sign = -1;
	}

	const double x = (M_PI / 2) * z * z;
	const double f = fresnel_f(z);
	const double g = fresnel_g(z);
	const double c = cos(x);
	const double s = sin(x);

	const double sum = 0.5 - (f * c) - (g * s);

	const double fc1 = (f + ERR_MARGIN) * c;
	const double fc2 = (f - ERR_MARGIN) * c;
	const double gs1 = (g + ERR_MARGIN) * s;
	const double gs2 = (g - ERR_MARGIN) * s;

	if (sign > 0)
	{
		low  =  (0.5 - max(fc1, fc2) - max(gs1, gs2));
		high =  (0.5 - min(fc1, fc2) - min(gs1, gs2));
	}
	else
	{
		low  = -(0.5 - min(fc1, fc2) - min(gs1, gs2));
		high = -(0.5 - max(fc1, fc2) - max(gs1, gs2));
	}

	return sign * sum;
}

double __fastcall E_Phase_Center(const double a, const double v, const double l, double &lv, double &hv)
{
	double fcl, fch;
	double fsl, fsh;

	const double q = a / (2 * v);

	const double fc = fresnel_cos(q, fcl, fch);
	const double fs = fresnel_sin(q, fsl, fsh);

	const double c2 = cos((M_PI / 2) * q * q);
	const double s2 = sin((M_PI / 2) * q * q);

	const double res = epc(l, q, c2, s2, fc, fs);

	hv = epc(l, q, c2, s2, fch, fsh);
	lv = epc(l, q, c2, s2, fcl, fsl);

	return res;
}

double __fastcall H_Phase_Center(const double a, const double v, const double l, double &lv, double &hv)
{
	double sul, suh;
	double swl, swh;
	double cul, cuh;
	double cwl, cwh;

	const double U  = (v / a) + (a / (2 * v));
	const double W  = (v / a) - (a / (2 * v));
	const double u2 = (M_PI / 2) * U * U;
	const double w2 = (M_PI / 2) * W * W;

	const double R = (W * cos(u2)) - (U * cos(w2));
	const double T = (U * sin(w2)) - (W * sin(u2));

	const double su = fresnel_sin(U, sul, suh);
	const double sw = fresnel_sin(W, swl, swh);
	const double cu = fresnel_cos(U, cul, cuh);
	const double cw = fresnel_cos(W, cwl, cwh);

	const double CD  = cu  - cw;
	const double CDL = cul - cwh;
	const double CDH = cuh - cwl;

	const double SD  = -su  + sw;
	const double SDL = -suh + swl;
	const double SDH = -sul + swh;

	const double res = hpc(l, R, T, CD,  SD);

	hv  = hpc(l, R, T, CDH, SDH);
	lv  = hpc(l, R, T, CDL, SDL);

	return res;
}

// convert to wavelengths for phase center routine
double __fastcall Phase_center(const char plane, const double aperture, const double guide, const double wavelength, const double Horn_length)
{
	double side_length;  // length of horn side in wavelengths
	double Ap_lambda;    // aperture in wavelengths
	double v;            // sqrt(side_length/2.0)
	double dmin, dmax;   // limits check on phase center
	double PC;           // phase center in wavelengths

	Ap_lambda   = aperture / wavelength;
	side_length = sqrt(Horn_length * Horn_length + ((aperture - guide) / 2) * ((aperture - guide) / 2));
	side_length = side_length / wavelength;
	v           = sqrt(side_length / 2);

	if (plane == 'E' || plane == 'e')
		PC = E_Phase_Center(Ap_lambda, v, side_length, dmin, dmax);
	else
		PC = H_Phase_Center(Ap_lambda, v, side_length, dmin, dmax);

	// test for bad phase center

	if (fabs(dmax - dmin) != (fabs(dmax) - fabs(dmin)))
	{
		//cerr << " \n\nERROR: BUG in " << plane << "_Phase_Center function\n\n";// call Matt\n";
	}
	else
	if (fabs(dmax - dmin) > (PC / 16) && (Ap_lambda > 1))
	{
		//cerr << " \n\nEXCESSIVE ERROR in " << plane << "_Phase_Center function \n\n";
		return 0;
	}

	return PC;
}

// compute the 'X' distance
double __fastcall computeX(const double f, const double y)
{
	return (y * y) / (4 * f);
}

// ************************************************************************
// dish beam width esitimation
//
// https://www.satsig.net/pointing/antenna-beamwidth-calculator.htm

// NB. this needs fixing to allow any efficiency value
double __fastcall factor_k(const double dia, const double Hz, const double eff)
{
	// based on edge illumination -3dB beamwidth    aperture eff   overall eff
	//          -infinity         72.8 lambda / D   0.75           0.6
	//          -10               65.3 lambda / D   0.83           0.674
	//          0                 58.4 lambda / D   1              0.8
	const double e8 = eff - 0.8;
	return 58.4 - (26.0 * e8) + (230.0 * e8 * e8);
}

double __fastcall half_power_beam_width(const double dia, const double Hz, const double k)
{
	const double lambda = sol / Hz;  // wave length
	return k * lambda / dia;
}

double __fastcall beam_width(const double half_power_beam_width, const double gain_reduction)
{
	return 2.0 * sqrt(half_power_beam_width * half_power_beam_width * gain_reduction * (1.0 / 12.0));
}

std::vector < std::pair <double, double> > __fastcall compute_beam_width(const double dia_meter, const double freq_Hz, const double eff_abs, const double max_dBi)
{
	const double k    = factor_k(dia_meter, freq_Hz, eff_abs);
	const double hpbw = half_power_beam_width(dia_meter, freq_Hz, k);

	std::vector < std::pair <double, double> > bw;

	if (dia_meter < 0.01 || freq_Hz < 1e5 || eff_abs < 0.01 || max_dBi < 1.0)
		return bw;

#if 0
	// table dB
	const double DB_TAB[] = {0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.2, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
	for (unsigned int i = 0; i < ARRAY_SIZE(DB_TAB); i++)
	{
		const double dB = DB_TAB[i];
		bw.push_back(std::make_pair(dB, beam_width(hpbw, dB)));
		// access the values ..
		//   bw[i].first;   // dB
		//   bw[i].second;  // beam width degress
		if (dB >= max_dBi)
			break;         // don't bother going behind the dish (backward radiation)
	}
	bw.push_back(std::make_pair(max_dBi, beam_width(hpbw, max_dBi)));
#else
	bool   _3dB = false;
	double dB   = 0.0;
	double step = max_dBi / 400;
	for (int i = 0; i < 100 && dB < max_dBi; i++)
	{
		if (!_3dB && dB >= 3.0)
		{	// include the -3dB beamwidth
				_3dB = true;
			dB = 3.0;
		}
		bw.push_back(std::make_pair(dB, beam_width(hpbw, dB)));
		dB += step;
		step *= 1.2;
	}
	bw.push_back(std::make_pair(max_dBi, beam_width(hpbw, max_dBi)));
#endif

	return bw;
}

// ************************************************************************

String __fastcall trimTrailingZeros(String s, const bool keep_dp)
{
	// remove trailing zeros
	while (!s.IsEmpty() && s[s.Length()] == '0')
		s = s.SubString(1, s.Length() - 1);

	// remove trailing decimal point
	if (!s.IsEmpty())
	{
		if (s[s.Length()] == '.')
		{
			if (!keep_dp)
				s = s.SubString(1, s.Length() - 1);
			else
				s += '0';
		}
	}

	return s;
}

bool __fastcall getBuildInfo(String filename, TVersion *version)
{
	DWORD ver_info_size;
	char *ver_info;
	UINT buffer_size;
	LPVOID buffer;
	DWORD dummy;

	if (version == NULL || filename.IsEmpty())
		return false;

	memset(version, 0, sizeof(TVersion));

	ver_info_size = ::GetFileVersionInfoSizeA(filename.c_str(), &dummy);
	if (ver_info_size == 0)
		return false;

	ver_info = new char [ver_info_size];
	if (ver_info == NULL)
		return false;

	if (::GetFileVersionInfoA(filename.c_str(), 0, ver_info_size, ver_info) == FALSE)
	{
		delete [] ver_info;
		return false;
	}

	if (::VerQueryValue(ver_info, _T("\\"), &buffer, &buffer_size) == FALSE)
	{
		delete [] ver_info;
		return false;
	}

	PVSFixedFileInfo ver = (PVSFixedFileInfo)buffer;
	version->MajorVer   = (ver->dwFileVersionMS >> 16) & 0xFFFF;
	version->MinorVer   = (ver->dwFileVersionMS >>  0) & 0xFFFF;
	version->ReleaseVer = (ver->dwFileVersionLS >> 16) & 0xFFFF;
	version->BuildVer   = (ver->dwFileVersionLS >>  0) & 0xFFFF;

	delete [] ver_info;

	return true;
}

// ************************************************************************

__fastcall TForm1::TForm1(TComponent* Owner)
	: TForm(Owner)
{
}

void __fastcall TForm1::FormCreate(TObject *Sender)
{
	{
//		char username[64];
//		DWORD size = sizeof(username);
//		if (::GetUserNameA(username, &size) != FALSE && size > 1)
//			m_ini_filename = ChangeFileExt(Application->ExeName, "_" + String(username) + ".ini");
//		else
			m_ini_filename = ChangeFileExt(Application->ExeName, ".ini");
	}

	::GetSystemInfo(&m_system_info);
//	sprintf(SystemInfoStr,
//					"OEM id: %u"crlf
//					"num of cpu's: %u"crlf
//					"page size: %u"crlf
//					"cpu type: %u"crlf
//					"min app addr: %lx"crlf
//					"max app addr: %lx"crlf
//					"active cpu mask: %u"crlf,
//					SystemInfo.dwOemId,
//					SystemInfo.dwNumberOfProcessors,
//					SystemInfo.dwPageSize,
//					SystemInfo.dwProcessorType,
//					SystemInfo.lpMinimumApplicationAddress,
//					SystemInfo.lpMaximumApplicationAddress,
//					SystemInfo.dwActiveProcessorMask);
//	MemoAddString(Memo1, SystemInfoStr);

	// the screen size is the phsyical screen size
	m_screen_width  = 0;
	m_screen_height = 0;
	HDC hDC = GetDC(0);
	if (hDC != NULL)
	{
		//ScreenBitsPerPixel = ::GetDeviceCaps(hDC, BITSPIXEL);
		m_screen_width  = ::GetDeviceCaps(hDC, HORZRES);
		m_screen_height = ::GetDeviceCaps(hDC, VERTRES);
		ReleaseDC(0, hDC);
	}

	{
		String s;
		TVersion version;
		getBuildInfo(Application->ExeName, &version);
		#ifdef _DEBUG
			s.printf("%s v%u.%u.%u.debug G6AMU", Application->Title.c_str(), version.MajorVer, version.MinorVer, version.ReleaseVer);
		#else
			s.printf("%s v%u.%u.%u G6AMU", Application->Title.c_str(), version.MajorVer, version.MinorVer, version.ReleaseVer);
		#endif
		this->Caption = s + "     compiled "__DATE__" "__TIME__;
	}

	this->DoubleBuffered   = true;
	Panel1->DoubleBuffered = true;
	Panel2->DoubleBuffered = true;
	Memo1->DoubleBuffered  = true;

	Memo1->Clear();
/*
	#ifdef __BORLANDC__
	{
		String s;
		s.printf("Borland %04x", __BORLANDC__);
		Memo1->Lines->Add(s);
	}
	#endif
*/
	m_mouse_x = -1;
	m_mouse_y = -1;

//	PaintBox1->Cursor = crNone;

	m_graph_bm = NULL;

	m_left_margin  = 0;
	m_right_margin = 0;

	m_pz_offset_x = 0.0f;
	m_pz_offset_y = 0.0f;
	m_pz_scale_x  = 1.0f;
	m_pz_scale_y  = 1.0f;

 	OpenDialog1->InitialDir = ExtractFilePath(Application->ExeName);
	SaveDialog1->InitialDir = ExtractFilePath(Application->ExeName);

	// ******************

	::PostMessage(this->Handle, WM_INIT_GUI, 0, 0);
}

void __fastcall TForm1::FormDestroy(TObject *Sender)
{
	//
}

void __fastcall TForm1::FormClose(TObject *Sender, TCloseAction &Action)
{
	saveSettings();

	if (m_graph_bm != NULL)
	{
		delete m_graph_bm;
		m_graph_bm = NULL;
	}
}

void __fastcall TForm1::WMWindowPosChanging(TWMWindowPosChanging &msg)
{
	const int thresh = 8;

	RECT work_area;
	SystemParametersInfo(SPI_GETWORKAREA, 0, &work_area, 0);

	const int dtLeft   = Screen->DesktopRect.left;
	const int dtRight  = Screen->DesktopRect.right;
	const int dtTop    = Screen->DesktopRect.top;
	const int dtBottom = Screen->DesktopRect.bottom;
	const int dtWidth  = dtRight - dtLeft;
	const int dtHeight = dtBottom - dtTop;

//	const int waLeft = work_area.left;
//	const int waTop = work_area.top;
//	const int waRight = work_area.right;
//	const int waBottom = work_area.bottom;
	const int waWidth = work_area.right - work_area.left;
	const int waHeight = work_area.bottom - work_area.top;

	int x = msg.WindowPos->x;
	int y = msg.WindowPos->y;
	int w = msg.WindowPos->cx;
	int h = msg.WindowPos->cy;

	{	// sticky screen edges
		if (std::abs((int)(x - work_area.left)) < thresh)
			x = work_area.left;			// stick left to left side
		else
		if (std::abs((int)((x + w) - work_area.right)) < thresh)
			x = work_area.right - w;	// stick right to right side

		if (std::abs((int)(y - work_area.top)) < thresh)
			y = work_area.top;			// stick top to top side
		else
		if (std::abs((int)((y + h) - work_area.bottom)) < thresh)
			y = work_area.bottom - h;	// stick bottom to bottm side

		// stick the right side to the right side of the screen if the left side is stuck to the left side of the screen
		if (x == work_area.left)
			if ((w >= (waWidth - thresh)) && (w <= (waWidth + thresh)))
				w = waWidth;

		// stick the bottom to the bottom of the screen if the top is stuck to the top of the screen
		if (y == work_area.top)
			if ((h >= (waHeight - thresh)) && (h <= (waHeight + thresh)))
				h = waHeight;
	}

	{	// limit minimum size
		if (w < Constraints->MinWidth)
			 w = Constraints->MinWidth;
		if (h < Constraints->MinHeight)
			 h = Constraints->MinHeight;
	}

	{	// limit maximum size
		if (w > Constraints->MaxWidth && Constraints->MaxWidth > Constraints->MinWidth)
			 w = Constraints->MaxWidth;
		if (h > Constraints->MaxHeight && Constraints->MaxHeight > Constraints->MinHeight)
			 h = Constraints->MaxHeight;
	}

	{	// limit maximum size
		if (w > dtWidth)
			 w = dtWidth;
		if (h > dtHeight)
			 h = dtHeight;
	}

	if (Application->MainForm && this != Application->MainForm)
	{	// stick to our main form sides
		const TRect rect = Application->MainForm->BoundsRect;

		if (std::abs((int)(x - rect.left)) < thresh)
			x = rect.left;			// stick to left to left side
		else
		if (std::abs((int)((x + w) - rect.left)) < thresh)
			x = rect.left - w;	// stick right to left side
		else
		if (std::abs((int)(x - rect.right)) < thresh)
			x = rect.right;		// stick to left to right side
		else
		if (std::abs((int)((x + w) - rect.right)) < thresh)
			x = rect.right - w;	// stick to right to right side

		if (std::abs((int)(y - rect.top)) < thresh)
			y = rect.top;			// stick top to top side
		else
		if (std::abs((int)((y + h) - rect.top)) < thresh)
			y = rect.top - h;		// stick bottom to top side
		else
		if (std::abs((int)(y - rect.bottom)) < thresh)
			y = rect.bottom;		// stick top to bottom side
		else
		if (std::abs((int)((y + h) - rect.bottom)) < thresh)
			y = rect.bottom - h;	// stick bottom to bottom side
	}

	{	// stop it completely leaving the desktop area
		if (x < (dtLeft - Width + (dtWidth / 15)))
			  x = dtLeft - Width + (dtWidth / 15);
		if (x > (dtWidth - (Screen->Width / 15)))
			  x = dtWidth - (Screen->Width / 15);
		if (y < dtTop)
			 y = dtTop;
		if (y > (dtBottom - (dtHeight / 10)))
			  y = dtBottom - (dtHeight / 10);
	}

	msg.WindowPos->x  = x;
	msg.WindowPos->y  = y;
	msg.WindowPos->cx = w;
	msg.WindowPos->cy = h;
}

void __fastcall TForm1::CMMouseEnter(TMessage &msg)
{
	TComponent *Comp = (TComponent *)msg.LParam;
	if (!Comp)
		return;

//	if (dynamic_cast<TControl *>(Comp) == NULL)
//		return;		// only interested in on screen controls

	if (dynamic_cast<TPaintBox *>(Comp) != NULL)
	{
		TPaintBox *pb = (TPaintBox *)Comp;
/*
		if (pb == PaintBox1)
		{
//			SetCursor(NULL);
			TCursor cursor = crNone;
			if (sender->Cursor != cursor)
			{
				pb->Cursor = cursor;
				pb->Parent->Perform(WM_SETCURSOR, (unsigned int)sender->Parent->Handle, MAKELPARAM(HTCLIENT, WM_MOUSEMOVE));
			}
		}
*/
		if (pb != PaintBox1 && (m_mouse_x >= 0 || m_mouse_y >= 0))
		{
			m_mouse_x = -1;
			m_mouse_y = -1;
			PaintBox1->Invalidate();
//			StatusBar1->Panels->Items[4]->Text = "";
//			StatusBar1->Update();
		}
	}
	else
	if (m_mouse_x >= 0 || m_mouse_y >= 0)
	{
		m_mouse_x = -1;
		m_mouse_y = -1;
		PaintBox1->Invalidate();
//		StatusBar1->Panels->Items[4]->Text = "";
//		StatusBar1->Update();
	}
}

void __fastcall TForm1::CMMouseLeave(TMessage &msg)
{
	TComponent *Comp = (TComponent *)msg.LParam;
	if (!Comp)
		return;

//	if (dynamic_cast<TControl *>(Comp) == NULL)
//		return;		// only interested in on screen controls

	if (dynamic_cast<TPaintBox *>(Comp) != NULL)
	{
		TPaintBox *pb = (TPaintBox *)Comp;
		if (pb == PaintBox1 || (m_mouse_x >= 0 || m_mouse_y >= 0))
		{
			m_mouse_x = -1;
			m_mouse_y = -1;
			PaintBox1->Invalidate();
//			StatusBar1->Panels->Items[4]->Text = "";
//			StatusBar1->Update();
		}
	}
	else
	if (m_mouse_x >= 0 || m_mouse_y >= 0)
	{
		m_mouse_x = -1;
		m_mouse_y = -1;
		PaintBox1->Invalidate();
//		StatusBar1->Panels->Items[4]->Text = "";
//		StatusBar1->Update();
	}
}

void __fastcall TForm1::WMInitGUI(TMessage &msg)
{
	loadSettings();

	BringToFront();
//	::SetForegroundWindow(Handle);

//	if (Application->MainForm)
//		Application->MainForm->Update();

	Memo1->Clear();

	::PostMessage(this->Handle, WM_COMPUTE_DISH, 0, 0);
}

#pragma option push
#pragma warn -8004   // silence "unused variable" warning

void __fastcall TForm1::WMComputeDish(TMessage &msg)
{
	//(void)msg;
	doUpdate();
}

#pragma option pop

void __fastcall TForm1::loadSettings()
{
	int          i;
	float        f;
	String       s;
	bool         b;
	TNotifyEvent ne;

	TIniFile *ini = new TIniFile(m_ini_filename);
	if (ini == NULL)
		return;

	Top    = ini->ReadInteger("MainForm", "Top",     Top);
	Left   = ini->ReadInteger("MainForm", "Left",   Left);
	Width  = ini->ReadInteger("MainForm", "Width",  Width);
	Height = ini->ReadInteger("MainForm", "Height", Height);

	FrequencyEdit->Text        = ini->ReadString( "params", "Frequency",  FrequencyEdit->Text.Trim());
#if 1
	DiameterEdit->Text         = ini->ReadString( "params", "Diameter",   DiameterEdit->Text.Trim());
#else
	GainEdit->Text             = ini->ReadString( "params", "Gain",       GainEdit->Text.Trim());
#endif
	FDRatioEdit->Text          = ini->ReadString( "params", "FDRatio",    FDRatioEdit->Text.Trim());
	EfficiencyCSpinEdit->Value = ini->ReadInteger("params", "Efficiency", EfficiencyCSpinEdit->Value);
	PointsCSpinEdit->Value     = ini->ReadInteger("params", "Points",     PointsCSpinEdit->Value);
	PetalsCSpinEdit->Value     = ini->ReadInteger("params", "Petals",     PetalsCSpinEdit->Value);

	delete ini;
}

void __fastcall TForm1::saveSettings()
{
	String s;

	DeleteFile(m_ini_filename);

	TIniFile *ini = new TIniFile(m_ini_filename);
	if (ini == NULL)
		return;

	ini->WriteInteger("MainForm", "Top",    Top);
	ini->WriteInteger("MainForm", "Left",   Left);
	ini->WriteInteger("MainForm", "Width",  Width);
	ini->WriteInteger("MainForm", "Height", Height);

	ini->WriteString( "params", "Frequency",  FrequencyEdit->Text.Trim());
#if 1
	ini->WriteString( "params", "Diameter",   DiameterEdit->Text.Trim());
#else
	ini->WriteString( "params", "Gain",       GainEdit->Text.Trim());
#endif
	ini->WriteString( "params", "FDRatio",    FDRatioEdit->Text.Trim());
	ini->WriteInteger("params", "Efficiency", EfficiencyCSpinEdit->Value);
	ini->WriteInteger("params", "Points",     PointsCSpinEdit->Value);
	ini->WriteInteger("params", "Petals",     PetalsCSpinEdit->Value);

	delete ini;
}

void __fastcall TForm1::comboBoxAutoWidth(TComboBox *comboBox)
{
	if (!comboBox)
		return;

	#define COMBOBOX_HORIZONTAL_PADDING	4

	int itemsFullWidth = comboBox->Width;

	// get the max needed with of the items in dropdown state
	for (int i = 0; i < comboBox->Items->Count; i++)
	{
		int itemWidth = comboBox->Canvas->TextWidth(comboBox->Items->Strings[i]);
		itemWidth += 2 * COMBOBOX_HORIZONTAL_PADDING;
		if (itemsFullWidth < itemWidth)
			itemsFullWidth = itemWidth;
	}

	if (comboBox->DropDownCount < comboBox->Items->Count)
		itemsFullWidth += ::GetSystemMetrics(SM_CXVSCROLL);

	::SendMessage(comboBox->Handle, CB_SETDROPPEDWIDTH, itemsFullWidth, 0);
}

void __fastcall TForm1::drawArc(Graphics::TBitmap *bm, const int x, const int y, const double start_deg, const double stop_deg, const double radius)
{
	const int x1 = x - radius;
	const int y1 = y - radius;
	const int x2 = x + radius;
	const int y2 = y + radius;

	const int x3 = x + (int)floor((sin((stop_deg) * deg_to_rad * 2) * radius) + 0.5);
	const int y3 = y - (int)floor((cos((stop_deg) * deg_to_rad * 2) * radius) + 0.5);

	const int x4 = x + (int)floor((sin((start_deg) * deg_to_rad * 2) * radius) + 0.5);
	const int y4 = y - (int)floor((cos((start_deg) * deg_to_rad * 2) * radius) + 0.5);

	bm->Canvas->Arc(x1, y1, x2, y2, x3, y3, x4, y4);
}

void __fastcall TForm1::drawCircle(Graphics::TBitmap *bm, const int x, const int y, const int radius)
{
	bm->Canvas->Ellipse(x - radius, y - radius, x + radius, y + radius);
}

void __fastcall TForm1::PaintBox1Paint(TObject *Sender)
{
	String s;

	TPaintBox *pb = dynamic_cast<TPaintBox *>(Sender);
	if (pb == NULL)
		return;

	if (m_graph_bm != NULL)
	{
		if (m_graph_bm->Width != pb->Width || m_graph_bm->Height != pb->Height)
		{
			delete m_graph_bm;
			m_graph_bm = NULL;
		}
	}

	if (m_graph_bm == NULL)
	{
		m_graph_bm = new Graphics::TBitmap();
		if (m_graph_bm == NULL)
		{
			pb->Canvas->Brush->Color = pb->Color;
			pb->Canvas->Brush->Style = bsSolid;
			pb->Canvas->FillRect(pb->Canvas->ClipRect);
			return;
		}

		m_graph_bm->Monochrome   = false;
		m_graph_bm->Transparent  = false;
		m_graph_bm->PixelFormat  = pf32bit;

		m_graph_bm->Canvas->Font = pb->Canvas->Font;
		m_graph_bm->Width        = pb->Width;
		m_graph_bm->Height       = pb->Height;
	}

	{
		m_graph_bm->Canvas->Pen->Style   = psSolid;
		m_graph_bm->Canvas->Pen->Color   = clGray;
		m_graph_bm->Canvas->Brush->Color = pb->Color;
		m_graph_bm->Canvas->Brush->Style = bsSolid;
		//m_graph_bm->Canvas->Rectangle(m_graph_bm->Canvas->ClipRect);
		m_graph_bm->Canvas->FillRect(m_graph_bm->Canvas->ClipRect);
	}

	const int th   = m_graph_bm->Canvas->TextHeight("!|qy");
	m_left_margin  = m_graph_bm->Canvas->TextWidth("  -00000  ");
	m_right_margin = m_graph_bm->Canvas->TextWidth("  -000  ");
	const int m_top_margin    = th * 2;
	const int m_bottom_margin = th * 2;

	const int w2 = m_graph_bm->Width  / 2;
	const int h2 = m_graph_bm->Height / 2;

	// **************

	if (m_frequency_Hz > 0 &&
		m_gain_dBi > 0.0 &&
		m_focus_diameter_ratio > 0.0 &&
		m_efficiency > 0.0 &&
		m_points > 0 &&
		!m_dish_point.empty())
	{
		const double dish_radius = m_dish_point[m_dish_point.size() - 1].r;
		const double dish_depth  = m_dish_point[m_dish_point.size() - 1].y;

		//const double front_to_focus = dish_depth - m_focus_point_meters;

		const int w = m_graph_bm->Width  - m_left_margin - m_right_margin;
		const int h = m_graph_bm->Height - m_top_margin  - m_bottom_margin;

		// ************************************
		// draw an estimated radiation pattern
		// doesn't include any side or rear lobes that would be present in the 'real' world

		if (!m_beam_width.empty())
		{
			double dB_max = m_beam_width[0].first;
			double dB_min = dB_max;

			for (unsigned int i = 1; i < m_beam_width.size(); i++)
			{
				const double dB = m_beam_width[i].first;
				if (dB_max < dB)
					dB_max = dB;
				if (dB_min > dB)
					dB_min = dB;
			}

			if (dB_min < dB_max && dB_max > 0.0)
			{
				const double scale = (double)w / dB_max;

				std::vector <TPoint> points1(m_beam_width.size());
				std::vector <TPoint> points2(m_beam_width.size());

				for (unsigned int i = 0; i < m_beam_width.size(); i++)
				{
					std::pair <double, double> bw = m_beam_width[i];

					const double dB  = bw.first;
					const double rad = (bw.second / 2) * deg_to_rad;

					const double sn  = sin(rad);
					const double cs  = cos(rad);

					const int x = (int)floor((cs * dB * scale) + 0.5);
					const int y = (int)floor((sn * dB * scale) + 0.5);

					points1[i] = TPoint(m_left_margin + x, h2 - y);
					points2[i] = TPoint(m_left_margin + x, h2 + y);
				}

				m_graph_bm->Canvas->Brush->Style = bsClear;
				//m_graph_bm->Canvas->Brush->Style = bsSolid;
				//m_graph_bm->Canvas->Brush->Color = pb->Color;
				m_graph_bm->Canvas->Pen->Color = (TColor)RGB(0,128, 0);
				//m_graph_bm->Canvas->Pen->Color = clGreen;
				m_graph_bm->Canvas->Pen->Style = psDot;
				m_graph_bm->Canvas->Pen->Width = 1;
				m_graph_bm->Canvas->Polyline(&points1[0], points1.size() - 1);
				m_graph_bm->Canvas->Polyline(&points2[0], points2.size() - 1);

				m_graph_bm->Canvas->Brush->Color = pb->Color;
				s.printf(" %0.1f dBi ", m_gain_dBi);
				const int x = (int)floor((m_gain_dBi * scale) + 0.5);
				m_graph_bm->Canvas->TextOut(m_left_margin + x - (m_graph_bm->Canvas->TextWidth(s) / 2), h2 - (th / 2), s);
			}
		}

		// ********************************************************
		// draw the petals/cutouts

		if (m_petals >= 3 && !m_petal_point.empty())
		{
			const double petal_radius = m_petal_point[m_petal_point.size() - 1].r;

			float scale = (h * 0.95f) / m_diameter_meters;
			if (scale > ((w * 0.95f) / (petal_radius * 2)))
				scale =   (w * 0.95f) / (petal_radius * 2);

			const int r = (int)floor((petal_radius * scale) + 0.5);

			const double petal_circumference = petal_radius * M_PI * 2;

			// draw a circle
			m_graph_bm->Canvas->Brush->Style = bsClear;
			//m_graph_bm->Canvas->Brush->Style = bsSolid;
			//m_graph_bm->Canvas->Brush->Color = pb->Color;
			m_graph_bm->Canvas->Pen->Color = clRed;
			m_graph_bm->Canvas->Pen->Style = psDot;
			m_graph_bm->Canvas->Pen->Width = 1;
			drawCircle(m_graph_bm, w2, h2, (int)floor((petal_radius * scale) + 0.5));
			{
				m_graph_bm->Canvas->Brush->Color = pb->Color;
				s.printf(" %0.3f cm radius ", petal_radius * 100);
				m_graph_bm->Canvas->TextOut(w2 - (m_graph_bm->Canvas->TextWidth(s) / 2), h2 - r - (th * 1), s);
				s.printf(" %0.3f cm circumference", petal_circumference * 100);
				m_graph_bm->Canvas->TextOut(w2 + 20, h2 - r + 5, s);
			}

			// draw the cutouts
			for (int pet = 0; pet < m_petals; pet++)
			{
				const double deg = (360.0 * pet) / m_petals;
				const double rad = deg * deg_to_rad;
				const double sn  = sin(rad);
				const double cs  = cos(rad);

				#if 0
				{	// petal/cutout center line
					m_graph_bm->Canvas->Brush->Style = bsClear;
					//m_graph_bm->Canvas->Brush->Style = bsSolid;
					//m_graph_bm->Canvas->Brush->Color = pb->Color;
					m_graph_bm->Canvas->Pen->Color = clDkGray;
					m_graph_bm->Canvas->Pen->Style = psDot;
					m_graph_bm->Canvas->Pen->Width = 1;

					const int x = (int)floor((sn * petal_radius * scale) + 0.5);
					const int y = (int)floor((cs * petal_radius * scale) + 0.5);

					m_graph_bm->Canvas->MoveTo(w2,     h2);
					m_graph_bm->Canvas->LineTo(w2 + x, h2 - y);
				}
				#endif

				m_graph_bm->Canvas->Brush->Style = bsClear;
				//m_graph_bm->Canvas->Brush->Style = bsSolid;
				//m_graph_bm->Canvas->Brush->Color = pb->Color;
				m_graph_bm->Canvas->Pen->Color = clRed;
				m_graph_bm->Canvas->Pen->Style = psSolid;
				m_graph_bm->Canvas->Pen->Width = 1;

				// right petal/cutout edge
				for (unsigned int i = 0; i < m_petal_point.size(); i++)
				{
					const t_petal_point p = m_petal_point[i];

					const double x = p.s * scale;
					//const double x = (p.w / 2) * scale;
					const double y = p.r * scale;

					const double xx = (cs * x) - (sn * y);
					const double yy = (sn * x) + (cs * y);

					const int xxx = (int)floor(xx + 0.5);
					const int yyy = (int)floor(yy + 0.5);

					if (i <= 0)
						m_graph_bm->Canvas->MoveTo(w2 + xxx, h2 - yyy);
					else
						m_graph_bm->Canvas->LineTo(w2 + xxx, h2 - yyy);
				}

				// left petal/cutout edge
				for (unsigned int i = 0; i < m_petal_point.size(); i++)
				{
					const t_petal_point p = m_petal_point[i];

					const double x = -p.s * scale;
					//const double x = -(p.w / 2) * scale;
					const double y = p.r * scale;

					const double xx = (cs * x) - (sn * y);
					const double yy = (sn * x) + (cs * y);

					const int xxx = (int)floor(xx + 0.5);
					const int yyy = (int)floor(yy + 0.5);

					if (i <= 0)
						m_graph_bm->Canvas->MoveTo(w2 + xxx, h2 - yyy);
					else
						m_graph_bm->Canvas->LineTo(w2 + xxx, h2 - yyy);
				}

				{	// draw the petal/cutout end edge
					const t_petal_point p = m_petal_point[m_petal_point.size() - 1];

					// distance from center
					const double y = p.r * scale;

					// left point
					const double x1 = -p.s * scale;
					//const double x1 = -(p.w / 2) * scale;

					// right point
					const double x2 =  p.s * scale;
					//const double x2 =  (p.w / 2) * scale;

					const double xx1 = (cs * x1) - (sn * y);
					const double yy1 = (sn * x1) + (cs * y);

					const double xx2 = (cs * x2) - (sn * y);
					const double yy2 = (sn * x2) + (cs * y);

					const int xxx1 = (int)floor(xx1 + 0.5);
					const int yyy1 = (int)floor(yy1 + 0.5);

					const int xxx2 = (int)floor(xx2 + 0.5);
					const int yyy2 = (int)floor(yy2 + 0.5);

					m_graph_bm->Canvas->MoveTo(w2 + xxx1, h2 - yyy1);
					m_graph_bm->Canvas->LineTo(w2 + xxx2, h2 - yyy2);
				}
			}
		}

		// ********************************************************
		// draw the dish and focus point

		{
			float scale = (h * 0.95f) / m_diameter_meters;
			if (scale > ((w * 0.95f) / m_focus_point_meters))
				scale =   (w * 0.95f) / m_focus_point_meters;
			if (scale > ((w * 0.95f) / m_dish_point[m_dish_point.size() - 1].y))
				scale =   (w * 0.95f) / m_dish_point[m_dish_point.size() - 1].y;

			const int dx = (int)floor((dish_depth  * scale) + 0.5);
			const int dy = (int)floor((dish_radius * scale) + 0.5);

			const int fx = (int)floor((m_focus_point_meters * scale) + 0.5);

			{	// v-line from dish edge to dish edge
				m_graph_bm->Canvas->Pen->Color = clBlack;
				m_graph_bm->Canvas->Pen->Style = psSolid;
				m_graph_bm->Canvas->Pen->Width = 1;
				m_graph_bm->Canvas->MoveTo(m_left_margin + dx, h2 - dy);
				m_graph_bm->Canvas->LineTo(m_left_margin + dx, h2 + dy);
			}

			{	// draw the dish curve

				m_graph_bm->Canvas->Brush->Style = bsClear;
				//m_graph_bm->Canvas->Brush->Style = bsSolid;
				//m_graph_bm->Canvas->Brush->Color = pb->Color;
				m_graph_bm->Canvas->Pen->Color = clBlack;
				m_graph_bm->Canvas->Pen->Style = psSolid;
				m_graph_bm->Canvas->Pen->Width = 1;

				int ox = 0;
				int oy = 0;

				for (unsigned int i = 0; i < m_dish_point.size(); i++)
				{
					const int x = (int)floor((m_dish_point[i].y * scale) + 0.5);
					const int y = (int)floor((m_dish_point[i].r * scale) + 0.5);

					// upper dot
					for (int y1 = -2; y1 <= 2; y1++)
						for (int x1 = -2; x1 <= 2; x1++)
							m_graph_bm->Canvas->Pixels[m_left_margin + x + x1][h2 - y + y1] = clBlack;

					// lower dot
					for (int y1 = -2; y1 <= 2; y1++)
						for (int x1 = -2; x1 <= 2; x1++)
							m_graph_bm->Canvas->Pixels[m_left_margin + x + x1][h2 + y + y1] = clBlack;

					// upper half
					m_graph_bm->Canvas->MoveTo(m_left_margin + ox, h2 - oy);
					m_graph_bm->Canvas->LineTo(m_left_margin + x,  h2 - y);

					// lower half
					m_graph_bm->Canvas->MoveTo(m_left_margin + ox, h2 + oy);
					m_graph_bm->Canvas->LineTo(m_left_margin + x,  h2 + y);

					ox = x;
					oy = y;
				}

				s.printf(" r %0.2f cm ", dish_radius * 100);
				m_graph_bm->Canvas->TextOut(m_left_margin + dx - (m_graph_bm->Canvas->TextWidth(s) / 2), h2 - dy - 3 - (th * 2), s);
				s.printf(" y %0.2f cm ", dish_depth * 100);
				m_graph_bm->Canvas->TextOut(m_left_margin + dx - (m_graph_bm->Canvas->TextWidth(s) / 2), h2 - dy - 3 - (th * 1), s);
			}

			{	// draw the horn
				const double horn_radius_meters = m_wave_length_meters * m_simple_horn_diameter / 2;

				m_graph_bm->Canvas->Brush->Style = bsClear;
				//m_graph_bm->Canvas->Brush->Style = bsSolid;
				//m_graph_bm->Canvas->Brush->Color = pb->Color;

				m_graph_bm->Canvas->Pen->Color = clDkGray;
				m_graph_bm->Canvas->Pen->Style = psDot;
				m_graph_bm->Canvas->Pen->Width = 2;

				const int y = (int)floor((horn_radius_meters * scale) + 0.5);

				m_graph_bm->Canvas->MoveTo(m_left_margin + fx - 30, h2 - y);
				m_graph_bm->Canvas->LineTo(m_left_margin + fx + 30, h2 - y);
				m_graph_bm->Canvas->MoveTo(m_left_margin + fx - 30, h2 + y);
				m_graph_bm->Canvas->LineTo(m_left_margin + fx + 30, h2 + y);
			}

			{	// draw the feed point

				const int cross_size = 7;

				//m_graph_bm->Canvas->Brush->Style = bsClear;
				//m_graph_bm->Canvas->Brush->Style = bsSolid;
				m_graph_bm->Canvas->Brush->Color = pb->Color;

				// show 2 lines at dish edges facing in direction of radiation
				m_graph_bm->Canvas->Pen->Color = clDkGray;
				m_graph_bm->Canvas->Pen->Style = psDot;
				m_graph_bm->Canvas->Pen->Width = 1;
				m_graph_bm->Canvas->MoveTo(m_left_margin + dx, h2 - dy);
				m_graph_bm->Canvas->LineTo(m_left_margin + fx, h2 - dy);
				m_graph_bm->Canvas->MoveTo(m_left_margin + dx, h2 + dy);
				m_graph_bm->Canvas->LineTo(m_left_margin + fx, h2 + dy);

				// v-line from dish edge to dish edge at focus point
				m_graph_bm->Canvas->Pen->Color = clDkGray;
				m_graph_bm->Canvas->Pen->Style = psDot;
				m_graph_bm->Canvas->Pen->Width = 1;
				m_graph_bm->Canvas->MoveTo(m_left_margin + fx, h2 - dy);
				m_graph_bm->Canvas->LineTo(m_left_margin + fx, h2 + dy);

				s.printf(" %0.2f cm ", (m_focus_point_meters - dish_depth) * 100);
				m_graph_bm->Canvas->TextOut(m_left_margin + (((dx + fx) - m_graph_bm->Canvas->TextWidth(s)) / 2), h2 + dy - th, s);

				// draw the feed angle
				drawArc(m_graph_bm, m_left_margin + fx, h2, 45 + ((m_feed_angle * rad_to_deg) / 2), 45 - ((m_feed_angle * rad_to_deg) / 2), 30);

				// show 2 lines from dish edge to focus point
				m_graph_bm->Canvas->Pen->Color = clDkGray;
				m_graph_bm->Canvas->Pen->Style = psDot;
				m_graph_bm->Canvas->Pen->Width = 1;
				m_graph_bm->Canvas->MoveTo(m_left_margin + fx, h2);
				m_graph_bm->Canvas->LineTo(m_left_margin + dx, h2 - dy);
				m_graph_bm->Canvas->MoveTo(m_left_margin + fx, h2);
				m_graph_bm->Canvas->LineTo(m_left_margin + dx, h2 + dy);

				// feed support arms
				s.printf(" %0.2f cm ", m_feed_arm_meters * 100);
				m_graph_bm->Canvas->TextOut(m_left_margin + (((dx + fx) - m_graph_bm->Canvas->TextWidth(s)) / 2), h2 - (dy / 2) - th, s);

				// h-line from dish center to focus
				m_graph_bm->Canvas->Pen->Color = clDkGray;
				m_graph_bm->Canvas->Pen->Style = psDot;
				m_graph_bm->Canvas->Pen->Width = 1;
				m_graph_bm->Canvas->MoveTo(m_left_margin,      h2);
				m_graph_bm->Canvas->LineTo(m_left_margin + fx, h2);

				s.printf(" %0.2f cm ", m_focus_point_meters * 100);
				m_graph_bm->Canvas->TextOut(m_left_margin + (((dx + fx) - m_graph_bm->Canvas->TextWidth(s)) / 2), h2 - th, s);

				// show a small cross at the feed point
				m_graph_bm->Canvas->Pen->Color = clBlack;
				m_graph_bm->Canvas->Pen->Style = psSolid;
				m_graph_bm->Canvas->Pen->Width = 2;
				m_graph_bm->Canvas->MoveTo(m_left_margin + fx - cross_size, h2);
				m_graph_bm->Canvas->LineTo(m_left_margin + fx + cross_size, h2);
				m_graph_bm->Canvas->MoveTo(m_left_margin + fx,              h2 - cross_size);
				m_graph_bm->Canvas->LineTo(m_left_margin + fx,              h2 + cross_size);

				s = "focus";
				m_graph_bm->Canvas->TextOut(m_left_margin + fx + 5, h2 - 2 - th, s);
				s.printf("%0.1f\xB0", m_feed_angle * rad_to_deg);
				m_graph_bm->Canvas->TextOut(m_left_margin + fx + 5, h2 + 2, s);
			}
		}

		// ************************************
	}

	// ********************************************************

	pb->Canvas->Lock();
	{
		pb->Canvas->CopyMode = cmSrcCopy;

		if (m_graph_bm != NULL)
			pb->Canvas->Draw(0, 0, m_graph_bm);

#if 0
		// draw the mouse point
		if (m_mouse_x >= m_left_margin &&
			m_mouse_x <= (pb->Width - m_right_margin) &&
			m_mouse_y >= m_top_margin &&
			m_mouse_y <= (pb->Height - m_bottom_margin))
		{
			//const int w = pb->Width - m_right_margin - m_left_margin;
			//const int h = pb->Height;

			const int x = m_mouse_x - m_left_margin;
			const int y = m_mouse_y;

			//const float xx = (float)x / w;
			//const float yy = 1.0f - ((float)y / h);

//			String s;
//			s.printf(" %0.1fHz  %-0.1f%s ", xx * (sample_rate_Hz / 2), min_y_scale + (yy * (max_y_scale - min_y_scale)), dB_scale ? "dB" : "");

			pb->Canvas->Pen->Style   = psDot;
			pb->Canvas->Pen->Color   = clBlack;
			pb->Canvas->Brush->Style = bsSolid;
			pb->Canvas->Brush->Color = pb->Color;
			pb->Canvas->MoveTo(m_left_margin + x, m_top_margin + 5);
			pb->Canvas->LineTo(m_left_margin + x, pb->Height - m_bottom_margin - 5);
			pb->Canvas->MoveTo(m_left_margin + 5, y);
			pb->Canvas->LineTo(pb->Width - m_right_margin - 5, y);
/*
			const int fs = pb->Canvas->Font->Size;
			pb->Canvas->Font->Size = fs + 3;
			const int tw = pb->Canvas->TextWidth(s);
			int tx = m_left_margin + x + 5;
			if (tx > (width - m_right_margin - tw))
				tx = width - m_right_margin - tw;
			pb->Canvas->Brush->Style = bsSolid;
			pb->Canvas->Brush->Color = pb->Color;
			pb->Canvas->TextOut(tx, y - 5 - pb->Canvas->TextHeight("|qg0"), s);
			pb->Canvas->Font->Size = fs;
*/
		}
#endif
	}

	pb->Canvas->Unlock();
}

void __fastcall TForm1::doUpdate()
{
	int    i;
	double d;
	String s;

	Memo1->Clear();

	m_dish_point.resize(0);
	m_petal_point.resize(0);
	m_beam_width.resize(0);

	// ***************
	// fetch user input

	m_frequency_Hz         = 0;
	m_diameter_meters      = 0;
	m_gain_dBi             = 0;
	m_focus_diameter_ratio = 0;
	m_efficiency           = 0;
	m_points               = 0;
	m_petals               = 0;

	if (!TryStrToFloat(FrequencyEdit->Text, d) || d < 0.001 || d > 500.0)
	{
		Memo1->Lines->Add("error: frequency (>= 0.001 <= 500)");
		return;
	}
	m_frequency_Hz = d * 1e9;

#if 1
	if (!TryStrToFloat(DiameterEdit->Text, d) || d <= 0.0 || d > 500000.0)
	{
		Memo1->Lines->Add("error: diameter (> 0.0 <= 500000)");
		return;
	}
	m_diameter_meters = d / 100;
#else
	if (!TryStrToFloat(GainEdit->Text, d) || d < 0.5 || d > 100.0)
	{
		Memo1->Lines->Add("error: gain (>= 0.5 <= 100)");
		return;
	}
	m_gain_dBi = d;
#endif

	if (!TryStrToFloat(FDRatioEdit->Text, d) || d < 0.18 || d > 5.0)
	{
		Memo1->Lines->Add("error: F/D ratio (>= 0.18 <= 5.0)");
		return;
	}
	m_focus_diameter_ratio = d;

	try
	{
		m_efficiency = (double)EfficiencyCSpinEdit->Value / 100;
	}
//	catch (Exception &exception)
//	{
//		Application->ShowException(&exception);
//	}
	catch (...)
	{
		Memo1->Lines->Add("error: efficiency (> 0 <= 100)");
		return;
	}

	try
	{
		m_points = PointsCSpinEdit->Value;
	}
	catch (...)
	{
		Memo1->Lines->Add("error: points (> 0)");
		return;
	}

	try
	{
		m_petals = PetalsCSpinEdit->Value;
	}
	catch (...)
	{
		Memo1->Lines->Add("error: petals (> 0)");
		return;
	}

	// *********************
	// compute the dish size/curve

	m_wave_length_meters = sol / m_frequency_Hz;

#if 1
	// gain to dish diameter
	m_gain = (m_diameter_meters * M_PI) / m_wave_length_meters;
	m_gain = (m_gain * m_gain) * m_efficiency;

	m_gain_dBi = 10.0f * log10(m_gain);

//	s.printf("%0.1f", m_gain_dBi);
//	GainEdit->Text = s;
#else
	// dish diameter to gain
	m_gain = pow(10.0, m_gain_dBi / 10.0);

	m_diameter_meters = (sqrt(m_gain / m_efficiency) * m_wave_length_meters) / M_PI;
//	m_diameter_meters = sqrt((m_gain * (m_wave_length_meters * m_wave_length_meters)) / (m_efficiency * (M_PI * M_PI)));

//	s.printf("%0.1f", m_diameter_meters * 100);
//	DiameterEdit->Text = s;
#endif

	m_focus_point_meters = m_diameter_meters * m_focus_diameter_ratio;

	for (int i = 0; i <= m_points; i++)
	{
		const double r = (m_diameter_meters * i) / (m_points * 2);
		const double y = (r * r) / (4 * m_focus_point_meters);
		m_dish_point.push_back(t_dish_point(r, y));
	}

	// *********************
	// compute the petal shape

	if (m_petals >= 3)
	{
		for (unsigned int i = 0; i < m_dish_point.size(); i++)
		{
			// https://www.instructables.com/How-to-build-a-strikeheliostatstrike-paraboli/

			const double n = m_petals;
			const double p = m_focus_point_meters;
			const double x = m_dish_point[i].r;

			const double y = (x * x) / (4 * p);
			const double r = (p * log(sqrt((x * x) + (4 * p * p)) + x)) + ((x / (4 * p)) * sqrt((x * x) + (4 * p * p))) - (p * log(2 * p));
			const double s = (M_PI / n) * (r - x);
			const double w = (2.0 * M_PI * r / n) - s;     // this is incorrect

			m_petal_point.push_back(t_petal_point(x, y, r, s, w));
		}
	}

	// *********************
	// compute feed antenna details

	m_feed_angle = 2.0 * atan2(1.0, (2.0 * m_focus_diameter_ratio) - (1.0 / (8.0 * m_focus_diameter_ratio)));

	m_feed_arm_meters = 2.0 * m_focus_point_meters / (1.0 + cos(m_feed_angle / 2));

	// space_attn is additional path loss to edge of dish
	m_space_attenuation_dB = (20.0 * log10(m_feed_arm_meters / m_focus_point_meters));

	// desired taper to achieve 10 dB actual edge taper
	m_desired_taper_dB = 10.0 - m_space_attenuation_dB;  // dB

	if (m_desired_taper_dB <= 0.0)
	{
		Memo1->Lines->Add("error: F/D ratio (>= 0.18)");
		return;
	}

	m_beam_width_3dB = m_feed_angle * sqrt(3.0 / m_desired_taper_dB);

	m_simple_horn_diameter = 66.0 / (m_beam_width_3dB * rad_to_deg);       // in wave length

	// *********************
	// estimate the dish beam width

	m_beam_width = compute_beam_width(m_diameter_meters, m_frequency_Hz, m_efficiency, m_gain_dBi);

	for (unsigned int i = 0; i < m_beam_width.size(); i++)
		m_beam_width[i].first = m_gain_dBi - m_beam_width[i].first;

	// *********************
	// show dish details

	Memo1->Lines->BeginUpdate();

	{
		Memo1->Lines->Add("");
		if (m_frequency_Hz < 1e3)
			s.printf("                 Frequency: %0.1g Hz", m_frequency_Hz / 1e0);
		else
		if (m_frequency_Hz < 1e6)
			s.printf("                 Frequency: %0.3g kHz", m_frequency_Hz / 1e3);
		else
		if (m_frequency_Hz < 1e9)
			s.printf("                 Frequency: %0.6g MHz", m_frequency_Hz / 1e6);
		else
		if (m_frequency_Hz < 1e12)
			s.printf("                 Frequency: %0.9g GHz", m_frequency_Hz / 1e9);
		else
			s.printf("                 Frequency: %0.12g THz", m_frequency_Hz / 1e12);
		Memo1->Lines->Add(s);
		s.printf("        Wave length/Lambda: %0.2f cm (%0.4f M)", m_wave_length_meters * 100, m_wave_length_meters);
		Memo1->Lines->Add(s);
		s.printf("                  Diameter: %0.1f cm (%0.3g M)", m_diameter_meters * 100, m_diameter_meters);
		Memo1->Lines->Add(s);
		s.printf("             Circumference: %0.1f cm (%0.3g M)", m_diameter_meters * M_PI * 100, m_diameter_meters * M_PI);
		Memo1->Lines->Add(s);
		s.printf("                      Gain: %0.2f dBi (%0.1f)", m_gain_dBi, m_gain);
		Memo1->Lines->Add(s);
		s.printf("Assumed overall efficiency: %0.1f%% (%0.2f)", m_efficiency * 100, m_efficiency);
		Memo1->Lines->Add(s);
		s.printf("          Feed focus point: %0.1f cm from back of dish", m_focus_point_meters * 100);
		Memo1->Lines->Add(s);
		s.printf("        Illumination angle: %0.2f\xB0", m_feed_angle * rad_to_deg);
		Memo1->Lines->Add(s);
		s.printf("         Space attenuation: %0.2f dB", m_space_attenuation_dB);
		Memo1->Lines->Add(s);
		s.printf("             Desired taper: %0.2f dB for 10dB edge illumination", m_desired_taper_dB);
		Memo1->Lines->Add(s);

		Memo1->Lines->Add("");
//		s.printf(" Simple feed horn diameter: %0.3f wavelengths ~ 3dB beamwidth of %0.2f\xB0", m_simple_horn_diameter, m_beam_width_3dB * rad_to_deg);
//		Memo1->Lines->Add(s);
		s.printf(" Simple feed horn diameter: %0.3f cm ~ 3dB beamwidth of %0.2f\xB0", m_wave_length_meters * m_simple_horn_diameter * 100, m_beam_width_3dB * rad_to_deg);
		Memo1->Lines->Add(s);
	}

	{
		Memo1->Lines->Add("");
		Memo1->Lines->Add("Dish curve");
		Memo1->Lines->Add("    x (radius)  y (from back)");
		for (unsigned int i = 0; i < m_dish_point.size(); i++)
		{
			const t_dish_point p = m_dish_point[i];
			s.printf(" %8.3fcm  %8.3fcm", p.r * 100, p.y * 100);
			Memo1->Lines->Add(s);
		}
	}

	if (m_petals >= 3 && !m_petal_point.empty())
	{
		Memo1->Lines->Add("");
		//Memo1->Lines->Add("Petal shape                                cut out     Petal");
		//Memo1->Lines->Add("    x (radius)  y (from back)  Radius      Segment     Width");
		Memo1->Lines->Add("Petal shape                                cut out");
		Memo1->Lines->Add("    x (radius)  y (from back)  Radius      Width");
		for (unsigned int i = 0; i < m_petal_point.size(); i++)
		{
			const t_petal_point p = m_petal_point[i];
			//s.printf(" %8.3fcm  %8.3fcm     %8.3fcm  %8.3fcm  %8.3fcm", p.x * 100, p.y * 100, p.r * 100, p.s * 100, p.w * 100);
			s.printf(" %8.3fcm  %8.3fcm     %8.3fcm  %8.3fcm", p.x * 100, p.y * 100, p.r * 100, p.s * 100);
			Memo1->Lines->Add(s);
		}
	}

	if (m_petals >= 3)
	{
		Memo1->Lines->Add("");
		Memo1->Lines->Add("Petal degrees");
		for (int i = 0; i < m_petals; )
		{
			s = "";
			for (int k = 0; k < 6 && i < m_petals; k++, i++)
			{
				String s2;
				s2.printf("%5.1f", (360.0 * i) / m_petals);
				s += s2 + "   ";
			}
			Memo1->Lines->Add(s.TrimRight());
		}
	}

	{	// show the beam width
		Memo1->Lines->Add("");
		Memo1->Lines->Add("Beam width");
		Memo1->Lines->Add("   dBi    beam width   rel dB");
		for (unsigned int i = 0; i < m_beam_width.size(); i++)
		{
			std::pair <double, double> bw = m_beam_width[i];
			const double dB  = bw.first;
			const double deg = bw.second;
			s.printf("  %6.2f   %5.1f       %6.2f", dB, deg, dB - m_gain_dBi);
			Memo1->Lines->Add(s);
		}
	}

	Memo1->Lines->EndUpdate();

	// ********************

	PaintBox1->Invalidate();
}

void __fastcall TForm1::FormKeyDown(TObject *Sender, WORD &Key,
		TShiftState Shift)
{
	switch (Key)
	{
		case VK_ESCAPE:
			Key = 0;
			Close();
			break;
		case VK_SPACE:
			Key = 0;
			break;
		case VK_UP:		// up arrow
//			Key = 0;
			break;
		case VK_DOWN:	// down arrow
//			Key = 0;
			break;
		case VK_LEFT:	// left arrow
//			Key = 0;
			break;
		case VK_RIGHT:	// right arrow
//			Key = 0;
			break;
		case VK_PRIOR:	// page up
//			Key = 0;
			break;
		case VK_NEXT:	// page down
//			Key = 0;
			break;
	}
}

// ****************************************************************************************

void __fastcall TForm1::PaintBox1MouseMove(TObject *Sender,
		TShiftState Shift, int X, int Y)
{
	m_mouse_x = X;
	m_mouse_y = Y;

//	String s;
//	s.printf("%-3d %-3d", m_mouse_x, m_mouse_y);
//	StatusBar1->Panels->Items[4]->Text = s;
//	StatusBar1->Update();

	PaintBox1->Invalidate();
}

void __fastcall TForm1::PaintBox1MouseDown(TObject *Sender,
		TMouseButton Button, TShiftState Shift, int X, int Y)
{
	//
}

void __fastcall TForm1::PaintBox1MouseUp(TObject *Sender,
		TMouseButton Button, TShiftState Shift, int X, int Y)
{
	//
}

void __fastcall TForm1::EditKeyDown(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
/*	switch (Key)
	{
		case VK_RETURN:
			Key = 0;
			doUpdate();
			break;
	}
*/
}

void __fastcall TForm1::EditKeyUp(TObject *Sender, WORD &Key,
		TShiftState Shift)
{
/*	switch (Key)
	{
		case VK_RETURN:
			Key = 0;
			doUpdate();
			break;
	}
*/
	doUpdate();
}

void __fastcall TForm1::CSpinChange(TObject *Sender)
{
	doUpdate();
}
