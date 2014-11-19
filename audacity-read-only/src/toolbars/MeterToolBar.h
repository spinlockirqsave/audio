/**********************************************************************

  Audacity: A Digital Audio Editor

  MeterToolbar.h

  Dominic Mazzoni
  Leland Lucius

  ToolBar to hold the VU Meter

**********************************************************************/

#ifndef __AUDACITY_METER_TOOLBAR__
#define __AUDACITY_METER_TOOLBAR__

#include "ToolBar.h"

class wxDC;
class wxGridBagSizer;
class wxSizeEvent;
class wxWindow;

class Meter;


// Constants used as bit pattern
const int kWithRecordMeter = 1;
const int kWithPlayMeter = 2;

class MeterToolBar:public ToolBar {

 public:

   MeterToolBar(int WhichMeters);
   virtual ~MeterToolBar();

   void Create(wxWindow *parent);
   bool DestroyChildren();

   virtual void Populate();
   virtual void Repaint(wxDC * WXUNUSED(dc)) {};
   virtual void EnableDisableButtons() {};
   virtual void UpdatePrefs();

   virtual void OnSize(wxSizeEvent & event);
   virtual bool Expose( bool show );

   int GetInitialWidth() {return (mWhichMeters == 
      (kWithRecordMeter + kWithPlayMeter)) ? 338 : 460;} // Separate bars used to be smaller.
   int GetMinToolbarWidth() { return 100; }
   wxSize GetDockedSize();

 private:
   void OnMeterPrefsUpdated(wxCommandEvent & evt);
   void RegenerateTooltips();

   int mWhichMeters;
   wxGridBagSizer *mSizer;
   Meter *mPlayMeter;
   Meter *mRecordMeter;

 public:

   DECLARE_CLASS(MeterToolBar);
   DECLARE_EVENT_TABLE();

};

namespace MeterToolBars {
   void AddMeters(Meter *playMeter, Meter *recordMeter);
   void RemoveMeters(Meter *playMeter, Meter *recordMeter);
   void GetMeters(Meter **playMeter, Meter **recordMeter);
   void StartMonitoring();
   void Clear();
   //Meter *mPlayMeter;
   //Meter *mRecordMeter;
};

#endif

