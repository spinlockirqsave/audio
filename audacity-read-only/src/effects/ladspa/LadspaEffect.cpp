/**********************************************************************

  Audacity: A Digital Audio Editor

  LadspaEffect.cpp

  Dominic Mazzoni

  This class implements a Ladspa Plug-in effect.

*******************************************************************//**

\class LadspaEffect
\brief An Effect that calls up a LADSPA plug in, i.e. many possible
effects from this one class.

*//****************************************************************//**

\class LadspaEffectDialog
\brief Dialog used with Effect

*//*******************************************************************/


#include "../../Audacity.h"

#include "ladspa.h"

#include <float.h>

#include <wx/wxprec.h>
#include <wx/button.h>
#include <wx/checkbox.h>
#include <wx/dynlib.h>
#include <wx/filename.h>
#include <wx/log.h>
#include <wx/menu.h>
#include <wx/msgdlg.h>
#include <wx/sizer.h>
#include <wx/slider.h>
#include <wx/statbox.h>
#include <wx/stattext.h>
#include <wx/textctrl.h>
#include <wx/tokenzr.h>
#include <wx/intl.h>
#include <wx/scrolwin.h>
#include <wx/version.h>

#include "LadspaEffect.h"       // This class's header file
#include "../../Internat.h"
#include "../../widgets/valnum.h"

// ============================================================================
// Module registration entry point
//
// This is the symbol that Audacity looks for when the module is built as a
// dynamic library.
//
// When the module is builtin to Audacity, we use the same function, but it is
// declared static so as not to clash with other builtin modules.
// ============================================================================
DECLARE_MODULE_ENTRY(AudacityModule)
{
   // Create and register the importer
   return new LadspaEffectsModule(moduleManager, path);
}

// ============================================================================
// Register this as a builtin module
// ============================================================================
DECLARE_BUILTIN_MODULE(LadspaBuiltin);

///////////////////////////////////////////////////////////////////////////////
//
// LadspaEffectsModule
//
///////////////////////////////////////////////////////////////////////////////

LadspaEffectsModule::LadspaEffectsModule(ModuleManagerInterface *moduleManager,
                                     const wxString *path)
{
   mModMan = moduleManager;
   if (path)
   {
      mPath = *path;
   }
}

LadspaEffectsModule::~LadspaEffectsModule()
{
}

// ============================================================================
// IdentInterface implementation
// ============================================================================

wxString LadspaEffectsModule::GetID()
{
   // Can be anything, but this is a v4 UUID
   return L"3ebd9fb9-c020-4c0d-a786-d6a914e55e31";
}

wxString LadspaEffectsModule::GetPath()
{
   return mPath;
}

wxString LadspaEffectsModule::GetName()
{
   return _("Ladspa Effect Module");
}

wxString LadspaEffectsModule::GetVendor()
{
   return _("The Audacity Team");
}

wxString LadspaEffectsModule::GetVersion()
{
   // This "may" be different if this were to be maintained as a separate DLL
   return LADSPAEFFECTS_VERSION;
}

wxString LadspaEffectsModule::GetDescription()
{
   return _("Provides Ladspa Effects");
}

// ============================================================================
// ModuleInterface implementation
// ============================================================================

bool LadspaEffectsModule::Initialize()
{
   // Nothing to do here
   return true;
}

void LadspaEffectsModule::Terminate()
{
   // Nothing to do here
   return;
}

bool LadspaEffectsModule::AutoRegisterPlugins(PluginManagerInterface & WXUNUSED(pm))
{
   return false;
}

wxArrayString LadspaEffectsModule::FindPlugins(PluginManagerInterface & pm)
{
#if defined(USE_LIBLRDF) && defined(EFFECT_CATEGORIES)

   EffectManager& em = EffectManager::Get();
   wxArrayString rdfPathList;
   wxString rdfPathVar;
   wxArrayString rdfFiles;

   InitCategoryMap();
   lrdf_init();

   rdfPathVar = wxGetenv(wxT("LADSPA_RDF_PATH"));
   if (rdfPathVar != wxT(""))
      wxGetApp().AddMultiPathsToPathList(rdfPathVar, rdfPathList);

#ifdef __WXGTK__
   wxGetApp().AddUniquePathToPathList(wxT("/usr/share/ladspa/rdf"),
                                      rdfPathList);
   wxGetApp().AddUniquePathToPathList(wxT("/usr/local/share/ladspa/rdf"),
                                      rdfPathList);
#endif

#ifdef __WXMAC__
   wxGetApp().AddUniquePathToPathList(wxT("/usr/share/ladspa/rdf"),
                                      rdfPathList);
   // XXX Maybe other Mac paths here?
#endif

#ifdef __WXMSW__
   //wxGetApp().AddUniquePathToPathList(wxT("WINDOWS LRDF PATH"),
   //                                   rdfPathList);
   // XXX Other Windows paths here.
#endif

   // Add the Audacity paths so we get ladspa.rdfs if we are using a local
   // liblrdf
   for(i=0; i<audacityPathList.GetCount(); i++) {
      wxString prefix = audacityPathList[i] + wxFILE_SEP_PATH;
      wxGetApp().AddUniquePathToPathList(prefix + wxT("rdf"),
                                         rdfPathList);
   }

   wxGetApp().FindFilesInPathList(wxT("*.rdf"), rdfPathList, rdfFiles);
   wxGetApp().FindFilesInPathList(wxT("*.rdfs"), rdfPathList, rdfFiles);
   for(size_t i = 0; i < rdfFiles.GetCount(); ++i) {
      wxString fileUri(wxT("file://"));
      fileUri += rdfFiles[i];
      lrdf_read_file(fileUri.mb_str(wxConvUTF8));
   }


   // Add all plugin categories found by LRDF
   lrdf_uris* cats =
      lrdf_get_all_subclasses("http://ladspa.org/ontology#Plugin");
   if (cats) {

      // Add the categories and find the plugins belonging to them
      for (size_t i = 0; i < cats->count; ++i) {
         char* label = lrdf_get_label(cats->items[i]);
         if (!label)
            continue;
         wxString uri = MapCategoryUri(wxString::FromAscii(cats->items[i]));
         em.AddCategory(uri, wxString::FromUTF8(label));
         std::free(label);

         lrdf_uris* plugs = lrdf_get_instances(cats->items[i]);
         if (plugs) {
            for (size_t j = 0; j < plugs->count; ++j) {
               unsigned long uid = lrdf_get_uid(plugs->items[j]);
               gPluginCategories.insert(std::make_pair(uid, uri));
            }
            lrdf_free_uris(plugs);
         }
      }

      // And their relationships
      for (size_t i = 0; i < cats->count; ++i) {
         EffectCategory* p =
            em.LookupCategory(MapCategoryUri(wxString::FromAscii(cats->
                                                                 items[i])));
         if (!p)
            continue;
         lrdf_uris* subs = lrdf_get_subclasses(cats->items[i]);
         if (subs) {
            for (size_t j = 0; j < subs->count; ++j) {
               EffectCategory* c =
                  em.LookupCategory(MapCategoryUri(wxString::FromAscii(subs->items[j])));
               if (c)
                  em.AddCategoryParent(c, p);
            }
            lrdf_free_uris(subs);
         }
      }

      lrdf_free_uris(cats);

   }

#endif

   wxArrayString pathList;
   wxArrayString files;
   wxString pathVar;

   // Check for the LADSPA_PATH environment variable
   pathVar = wxString::FromUTF8(getenv("LADSPA_PATH"));
   if (!pathVar.empty())
   {
      wxStringTokenizer tok(pathVar);
      while (tok.HasMoreTokens())
      {
         pathList.Add(tok.GetNextToken());
      }
   }

#if defined(__WXMAC__)
#define LADSPAPATH wxT("/Library/Audio/Plug-Ins/LADSPA")

   // Look in ~/Library/Audio/Plug-Ins/LADSPA and /Library/Audio/Plug-Ins/LADSPA
   pathList.Add(wxGetHomeDir() + wxFILE_SEP_PATH + LADSPAPATH);
   pathList.Add(LADSPAPATH);

   // Recursively scan for all shared objects
   pm.FindFilesInPathList(wxT("*.so"), pathList, files, true);

#elif defined(__WXMSW__)

   // Recursively scan for all DLLs
   pm.FindFilesInPathList(wxT("*.dll"), pathList, files, true);

#else
   
   pathList.Add(wxGetHomeDir() + wxFILE_SEP_PATH + wxT(".ladspa"));
   pathList.Add(wxT("/usr/local/lib/ladspa"));
   pathList.Add(wxT("/usr/lib/ladspa"));
   pathList.Add(wxT(LIBDIR) wxT("/ladspa"));

   // Recursively scan for all shared objects
   pm.FindFilesInPathList(wxT("*.so"), pathList, files, true);

#endif

   return files;
}

bool LadspaEffectsModule::RegisterPlugin(PluginManagerInterface & pm, const wxString & path)
{
   // Since we now have builtin VST support, ignore the VST bridge as it
   // causes duplicate menu entries to appear.
   wxFileName f(path);
   if (f.GetName().CmpNoCase(wxT("vst-bridge")) == 0) {
      return false;
   }

   // As a courtesy to some plug-ins that might be bridges to
   // open other plug-ins, we set the current working
   // directory to be the plug-in's directory.

   wxString saveOldCWD = ::wxGetCwd();
   wxString prefix = ::wxPathOnly(path);
   ::wxSetWorkingDirectory(prefix);

   int index = 0;
   LADSPA_Descriptor_Function mainFn = NULL;
   wxDynamicLibrary lib;
   if (lib.Load(path, wxDL_NOW)) {
      wxLogNull logNo;

      mainFn = (LADSPA_Descriptor_Function) lib.GetSymbol(wxT("ladspa_descriptor"));
      if (mainFn) {
         const LADSPA_Descriptor *data;

         for (data = mainFn(index); data; data = mainFn(++index)) {
#if defined(USE_LIBLRDF) && defined(EFFECT_CATEGORIES)
            std::set<wxString> categories;
            std::multimap<unsigned long, wxString>::const_iterator iter;
            iter = gPluginCategories.lower_bound(data->UniqueID);
            for ( ; (iter != gPluginCategories.end() &&
                     iter->first == data->UniqueID); ++iter)
               categories.insert(iter->second);
#endif
            LadspaEffect effect(path, index);
            if (effect.SetHost(NULL)) {
               pm.RegisterEffectPlugin(this, &effect);
            }
         }
      }
   }

   if (lib.IsLoaded()) {
      lib.Unload();
   }

   ::wxSetWorkingDirectory(saveOldCWD);

   return index > 0;
}

bool LadspaEffectsModule::IsPluginValid(const PluginID & ID,
                                        const wxString & path)
{
   return wxFileName::FileExists(path);
}

IdentInterface *LadspaEffectsModule::CreateInstance(const PluginID & ID,
                                                    const wxString & path)
{
   // For us, the ID is two words.
   // 1)  The Ladspa descriptor index
   // 2)  The library's path
   long index;
   ID.BeforeFirst(wxT(' ')).ToLong(&index);

   return new LadspaEffect(path, (int) index);
}

void LadspaEffectsModule::DeleteInstance(IdentInterface *instance)
{
   LadspaEffect *effect = dynamic_cast<LadspaEffect *>(instance);
   if (effect)
   {
      delete effect;
   }
}

///////////////////////////////////////////////////////////////////////////////
//
// LadspaEffectEventHelper
//
///////////////////////////////////////////////////////////////////////////////

enum
{
   ID_DURATION = 20000,
   ID_TOGGLES = 21000,
   ID_SLIDERS = 22000,
   ID_TEXTS = 23000,
};

BEGIN_EVENT_TABLE(LadspaEffectEventHelper, wxEvtHandler)
   EVT_COMMAND_RANGE(ID_TOGGLES, ID_TOGGLES + 999, wxEVT_COMMAND_CHECKBOX_CLICKED, LadspaEffectEventHelper::OnCheckBox)
   EVT_COMMAND_RANGE(ID_SLIDERS, ID_SLIDERS + 999, wxEVT_COMMAND_SLIDER_UPDATED, LadspaEffectEventHelper::OnSlider)
   EVT_COMMAND_RANGE(ID_TEXTS, ID_TEXTS + 999, wxEVT_COMMAND_TEXT_UPDATED, LadspaEffectEventHelper::OnTextCtrl)
END_EVENT_TABLE()

LadspaEffectEventHelper::LadspaEffectEventHelper(LadspaEffect *effect)
{
   mEffect = effect;
}

LadspaEffectEventHelper::~LadspaEffectEventHelper()
{
}

// ============================================================================
// LadspaEffectEventHelper implementation
// ============================================================================

void LadspaEffectEventHelper::OnCheckBox(wxCommandEvent & evt)
{
   mEffect->OnCheckBox(evt);
}

void LadspaEffectEventHelper::OnSlider(wxCommandEvent & evt)
{
   mEffect->OnSlider(evt);
}

void LadspaEffectEventHelper::OnTextCtrl(wxCommandEvent & evt)
{
   mEffect->OnTextCtrl(evt);
}

///////////////////////////////////////////////////////////////////////////////
//
// LadspaEffect
//
///////////////////////////////////////////////////////////////////////////////

LadspaEffect::LadspaEffect(const wxString & path, int index)
{
   mPath = path;
   mIndex = index;
   mData = NULL;

   mHost = NULL;
   mMaster = NULL;
   mReady = false;

   mInteractive = false;

   mAudioIns = 0;
   mAudioOuts = 0;
   mNumInputControls = 0;
   mNumOutputControls = 0;
   mSampleRate = 44100;

   mInputPorts = NULL;
   mOutputPorts = NULL;
   mInputControls = NULL;
   mOutputControls = NULL;

   mLatencyPort = -1;

   mEventHelper = NULL;
   mDialog = NULL;
   mSliders = NULL;
   mFields = NULL;
   mLabels = NULL;
   mToggles = NULL;
}

LadspaEffect::~LadspaEffect()
{
   if (mInputPorts)
   {
      delete [] mInputPorts;
   }

   if (mOutputPorts)
   {
      delete [] mOutputPorts;
   }

   if (mInputControls)
   {
      delete [] mInputControls;
   }

   if (mOutputControls)
   {
      delete [] mOutputControls;
   }

   if (mToggles)
   {
      delete [] mToggles;
   }

   if (mSliders)
   {
      delete [] mSliders;
   }

   if (mFields)
   {
      delete [] mFields;
   }

   if (mLabels)
   {
      delete [] mLabels;
   }
}

// ============================================================================
// IdentInterface implementation
// ============================================================================

wxString LadspaEffect::GetID()
{
   return wxString::Format(wxT("%d %s"), mIndex, mPath.c_str());
}

wxString LadspaEffect::GetPath()
{
   return mPath;
}

wxString LadspaEffect::GetName()
{
   return LAT1CTOWX(mData->Name);
}

wxString LadspaEffect::GetVendor()
{
   return LAT1CTOWX(mData->Maker);
}

wxString LadspaEffect::GetVersion()
{
   return _("N/A");
}

wxString LadspaEffect::GetDescription()
{
   return LAT1CTOWX(mData->Copyright);
}

// ============================================================================
// EffectIdentInterface implementation
// ============================================================================

EffectType LadspaEffect::GetType()
{
   if (mAudioIns == 0 && mAudioOuts == 0)
   {
      return EffectTypeNone;
   }

   if (mAudioIns == 0)
   {
      return EffectTypeGenerate;
   }

   if (mAudioOuts == 0)
   {
      return EffectTypeAnalyze;
   }

   return EffectTypeProcess;
}

wxString LadspaEffect::GetFamily()
{
   return LADSPAEFFECTS_FAMILY;
}

bool LadspaEffect::IsInteractive()
{
   return mInteractive;
}

bool LadspaEffect::IsDefault()
{
   return false;
}

bool LadspaEffect::IsLegacy()
{
   return false;
}

bool LadspaEffect::SupportsRealtime()
{
   return GetType() == EffectTypeProcess;
}

bool LadspaEffect::SupportsAutomation()
{
   return mNumInputControls > 0;
}

// ============================================================================
// EffectClientInterface Implementation
// ============================================================================

bool LadspaEffect::SetHost(EffectHostInterface *host)
{
   mHost = host;

   if (!Load())
   {
      return false;
   }

   mInputPorts = new unsigned long [mData->PortCount];
   mOutputPorts = new unsigned long [mData->PortCount];
   mInputControls = new float [mData->PortCount];
   mOutputControls = new float [mData->PortCount];

   for (unsigned long p = 0; p < mData->PortCount; p++)
   {
      LADSPA_PortDescriptor d = mData->PortDescriptors[p];

      // Collect the audio ports
      if (LADSPA_IS_PORT_AUDIO(d))
      {
         if (LADSPA_IS_PORT_INPUT(d)) 
         {
            mInputPorts[mAudioIns++] = p;
         }
         else if (LADSPA_IS_PORT_OUTPUT(d))
         {
            mOutputPorts[mAudioOuts++] = p;
         }
      }
      // Determine the port's default value
      else if (LADSPA_IS_PORT_CONTROL(d) && LADSPA_IS_PORT_INPUT(d))
      {
         mInteractive = true;

         LADSPA_PortRangeHint hint = mData->PortRangeHints[p];
         float val = float(1.0);
         float lower = hint.LowerBound;
         float upper = hint.UpperBound;

         if (LADSPA_IS_HINT_SAMPLE_RATE(hint.HintDescriptor))
         {
            lower *= mSampleRate;
            upper *= mSampleRate;
         }

         if (LADSPA_IS_HINT_BOUNDED_BELOW(hint.HintDescriptor) && val < lower)
         {
            val = lower;
         }

         if (LADSPA_IS_HINT_BOUNDED_ABOVE(hint.HintDescriptor) && val > upper)
         {
            val = upper;
         }

         if (LADSPA_IS_HINT_DEFAULT_MINIMUM(hint.HintDescriptor))
         {
            val = lower;
         }

         if (LADSPA_IS_HINT_DEFAULT_MAXIMUM(hint.HintDescriptor))
         {
            val = upper;
         }

         if (LADSPA_IS_HINT_DEFAULT_LOW(hint.HintDescriptor))
         {
            if (LADSPA_IS_HINT_LOGARITHMIC(hint.HintDescriptor))
            {
               val = exp(log(lower)) * 0.75f + log(upper) * 0.25f;
            }
            else
            {
               val = lower * 0.75f + upper * 0.25f;
            }
         }

         if (LADSPA_IS_HINT_DEFAULT_MIDDLE(hint.HintDescriptor))
         {
            if (LADSPA_IS_HINT_LOGARITHMIC(hint.HintDescriptor))
            {
               val = exp(log(lower)) * 0.5f + log(upper) * 0.5f;
            }
            else
            {
               val = lower * 0.5f + upper * 0.5f;
            }
         }

         if (LADSPA_IS_HINT_DEFAULT_HIGH(hint.HintDescriptor))
         {
            if (LADSPA_IS_HINT_LOGARITHMIC(hint.HintDescriptor))
            {
               val = exp(log(lower)) * 0.25f + log(upper) * 0.75f;
            }
            else
            {
               val = lower * 0.25f + upper * 0.75f;
            }
         }

         if (LADSPA_IS_HINT_DEFAULT_0(hint.HintDescriptor))
         {
            val = 0.0f;
         }

         if (LADSPA_IS_HINT_DEFAULT_1(hint.HintDescriptor))
         {
            val = 1.0f;
         }

         if (LADSPA_IS_HINT_DEFAULT_100(hint.HintDescriptor))
         {
            val = 100.0f;
         }

         if (LADSPA_IS_HINT_DEFAULT_440(hint.HintDescriptor))
         {
            val = 440.0f;
         }

         mNumInputControls++;
         mInputControls[p] = val;
      }
      else if (LADSPA_IS_PORT_CONTROL(d) && LADSPA_IS_PORT_OUTPUT(d))
      {
         mInteractive = true;

         mNumOutputControls++;
         mOutputControls[p] = 0.0;
 
         // Ladspa effects have a convention of providing latency on an output
         // control port whose name is "latency".
         if (strcmp(mData->PortNames[p], "latency") == 0)
         {
            mLatencyPort = p;
         }
      }
   }

   // mHost will be null during registration
   if (mHost)
   {
      mHost->GetSharedConfig(wxT("Settings"), wxT("BufferSize"), mUserBlockSize, 8192);
      mBlockSize = mUserBlockSize;

      bool haveDefaults;
      mHost->GetPrivateConfig(wxT("Default"), wxT("Initialized"), haveDefaults, false);
      if (!haveDefaults)
      {
         SaveParameters(wxT("Default"));
         mHost->SetPrivateConfig(wxT("Default"), wxT("Initialized"), true);
      }

      LoadParameters(wxT("Current"));
   }

   return true;
}

int LadspaEffect::GetAudioInCount()
{
   return mAudioIns;
}

int LadspaEffect::GetAudioOutCount()
{
   return mAudioOuts;
}

int LadspaEffect::GetMidiInCount()
{
   return 0;
}

int LadspaEffect::GetMidiOutCount()
{
   return 0;
}

void LadspaEffect::SetSampleRate(sampleCount rate)
{
   mSampleRate = rate;
}

sampleCount LadspaEffect::GetBlockSize(sampleCount maxBlockSize)
{
#if 0
   // TODO:  Allow user to specify the max blocksize
   if (mUserBlockSize > maxBlockSize)
   {
      mBlockSize = maxBlockSize;
   }
   else
   {
      mBlockSize = mUserBlockSize;
   }
#endif

   return mBlockSize;
}

sampleCount LadspaEffect::GetLatency()
{
   if (mLatencyPort >= 0 && !mLatencyDone)
   {
      mLatencyDone = true;
      return mOutputControls[mLatencyPort] * 2;
   }

   return 0;
}

sampleCount LadspaEffect::GetTailSize()
{
   return 0;
}

bool LadspaEffect::IsReady()
{
   return mReady;
}

bool LadspaEffect::ProcessInitialize()
{
   /* Instantiate the plugin */
   if (!mReady)
   {
      mMaster = InitInstance(mSampleRate);
      if (!mMaster)
      {
         return false;
      }
      mReady = true;
   }

   mLatencyDone = false;

   return true;
}

bool LadspaEffect::ProcessFinalize()
{
   if (mReady)
   {
      mReady = false;

      FreeInstance(mMaster);
      mMaster = NULL;
   }

   return true;
}

sampleCount LadspaEffect::ProcessBlock(float **inbuf, float **outbuf, sampleCount size)
{
   for (int i = 0; i < mAudioIns; i++)
   {
      mData->connect_port(mMaster, mInputPorts[i], inbuf[i]);
   }

   for (int i = 0; i < mAudioOuts; i++)
   {
      mData->connect_port(mMaster, mOutputPorts[i], outbuf[i]);
   }

   mData->run(mMaster, size);

   RefreshControls(true);

   return size;
}

bool LadspaEffect::RealtimeInitialize()
{
   return true;
}

bool LadspaEffect::RealtimeFinalize()
{
   for (size_t i = 0, cnt = mSlaves.GetCount(); i < cnt; i++)
   {
      FreeInstance(mSlaves[i]);
   }
   mSlaves.Clear();

   return true;
}

bool LadspaEffect::RealtimeSuspend()
{
   return true;
}

bool LadspaEffect::RealtimeResume()
{
   return true;
}

sampleCount LadspaEffect::RealtimeProcess(int group,
                                          float **inbuf,
                                          float **outbuf,
                                          sampleCount numSamples)
{
   if (group < 0 || group >= (int) mSlaves.GetCount())
   {
      return 0;
   }

   for (int i = 0; i < mAudioIns; i++)
   {
      mData->connect_port(mSlaves[group], mInputPorts[i], inbuf[i]);
   }

   for (int i = 0; i < mAudioOuts; i++)
   {
      mData->connect_port(mSlaves[group], mOutputPorts[i], outbuf[i]);
   }

   mData->run(mSlaves[group], numSamples);

   return numSamples;
}

bool LadspaEffect::RealtimeAddProcessor(int numChannels, float sampleRate)
{
   LADSPA_Handle slave = InitInstance(sampleRate);
   if (!slave)
   {
      return false;
   }

   mSlaves.Add(slave);
   mSlaveChannels.Add(numChannels);

   return true;
}

bool LadspaEffect::ShowInterface(wxWindow *parent, bool forceModal)
{
   if (mDialog)
   {
      mDialog->Close(true);
      return false;
   }

   mDialog = mHost->CreateUI(parent, this);
   if (!mDialog)
   {
      return false;
   }

   if ((SupportsRealtime() || GetType() == EffectTypeAnalyze) && !forceModal)
   {
      mDialog->Show();

      return false;
   }

   bool res = mDialog->ShowModal() != 0;
   mDialog = NULL;

   return res;
}

bool LadspaEffect::GetAutomationParameters(EffectAutomationParameters & parms)
{
   for (unsigned long p = 0; p < mData->PortCount; p++)
   {
      LADSPA_PortDescriptor d = mData->PortDescriptors[p];

      if (LADSPA_IS_PORT_CONTROL(d) && LADSPA_IS_PORT_INPUT(d))
      {
         if (!parms.Write(LAT1CTOWX(mData->PortNames[p]), mInputControls[p]))
         {
            return false;
         }
      }
   }

   return true;
}

bool LadspaEffect::SetAutomationParameters(EffectAutomationParameters & parms)
{
   for (unsigned long p = 0; p < mData->PortCount; p++)
   {
      LADSPA_PortDescriptor d = mData->PortDescriptors[p];

      if (LADSPA_IS_PORT_CONTROL(d) && LADSPA_IS_PORT_INPUT(d))
      {
         wxString labelText = LAT1CTOWX(mData->PortNames[p]);
         double d = 0.0;
         if (!parms.Read(labelText, &d))
         {
            return false;
         }

         mInputControls[p] = d;
      }
   }

   return true;
}

// ============================================================================
// EffectUIClientInterface Implementation
// ============================================================================

void LadspaEffect::SetUIHost(EffectUIHostInterface *host)
{
   mUIHost = host;
}

bool LadspaEffect::PopulateUI(wxWindow *parent)
{
   mParent = parent;

   mEventHelper = new LadspaEffectEventHelper(this);
   mParent->PushEventHandler(mEventHelper);

   mToggles = new wxCheckBox*[mData->PortCount];
   mSliders = new wxSlider*[mData->PortCount];
   mFields = new wxTextCtrl*[mData->PortCount];
   mLabels = new wxStaticText*[mData->PortCount];

   memset(mFields, 0, mData->PortCount * sizeof(wxTextCtrl *));

   wxSizer *marginSizer = new wxBoxSizer(wxVERTICAL);

   if (mNumInputControls)
   {
      wxSizer *paramSizer = new wxStaticBoxSizer(wxVERTICAL, mParent, _("Effect Settings"));

      wxFlexGridSizer *gridSizer = new wxFlexGridSizer(5, 0, 0);
      gridSizer->AddGrowableCol(3);

      wxControl *item;

      // Add the duration control for generators
      if (GetType() == EffectTypeGenerate)
      {
         item = new wxStaticText(mParent, 0, _("Duration:"));
         gridSizer->Add(item, 0, wxALIGN_CENTER_VERTICAL | wxALIGN_RIGHT | wxALL, 5);
         mDuration = new NumericTextCtrl(NumericConverter::TIME,
                                        mParent,
                                        ID_DURATION,
                                        _("hh:mm:ss + milliseconds"),
                                        mHost->GetDuration(),
                                        mSampleRate,
                                        wxDefaultPosition,
                                        wxDefaultSize,
                                        true);
         mDuration->SetName(_("Duration"));
         mDuration->EnableMenu();
         gridSizer->Add(mDuration, 0, wxALIGN_CENTER_VERTICAL | wxALL, 5);
         gridSizer->Add(1, 1, 0);
         gridSizer->Add(1, 1, 0);
         gridSizer->Add(1, 1, 0);
      }

      for (unsigned long p = 0; p < mData->PortCount; p++)
      {
         LADSPA_PortDescriptor d = mData->PortDescriptors[p];
         if (LADSPA_IS_PORT_AUDIO(d) || LADSPA_IS_PORT_OUTPUT(d))
         {
            continue;
         }

         wxString labelText = LAT1CTOWX(mData->PortNames[p]);
         item = new wxStaticText(mParent, 0, labelText + wxT(":"));
         gridSizer->Add(item, 0, wxALIGN_CENTER_VERTICAL | wxALIGN_RIGHT | wxALL, 5);

         wxString fieldText;
         LADSPA_PortRangeHint hint = mData->PortRangeHints[p];

         if (LADSPA_IS_HINT_TOGGLED(hint.HintDescriptor))
         {
            mToggles[p] = new wxCheckBox(mParent, ID_TOGGLES + p, wxT(""));
            mToggles[p]->SetName(labelText);
            mToggles[p]->SetValue(mInputControls[p] > 0);
            gridSizer->Add(mToggles[p], 0, wxALL, 5);

            gridSizer->Add(1, 1, 0);
            gridSizer->Add(1, 1, 0);
            gridSizer->Add(1, 1, 0);
            continue;
         }

         wxString bound;
         double lower = -FLT_MAX;
         double upper = FLT_MAX;
         bool haslo = false;
         bool hashi = false;
         bool forceint = false;

         if (LADSPA_IS_HINT_BOUNDED_BELOW(hint.HintDescriptor))
         {
            lower = hint.LowerBound;
            haslo = true;
         }

         if (LADSPA_IS_HINT_BOUNDED_ABOVE(hint.HintDescriptor))
         {
            upper = hint.UpperBound;
            hashi = true;
         }

         if (LADSPA_IS_HINT_SAMPLE_RATE(hint.HintDescriptor))
         {
            lower *= mSampleRate;
            upper *= mSampleRate;
            forceint = true;
         }

         // Don't specify a value at creation time.  This prevents unwanted events
         // being sent to the OnTextCtrl() handler before the associated slider
         // has been created.
         mFields[p] = new wxTextCtrl(mParent, ID_TEXTS + p);
         mFields[p]->SetName(labelText);
         gridSizer->Add(mFields[p], 0, wxALIGN_CENTER_VERTICAL | wxALL, 5);

         wxString str;
         if (haslo)
         {
            if (LADSPA_IS_HINT_INTEGER(hint.HintDescriptor) || forceint)
            {
               str.Printf(wxT("%d"), (int)(lower + 0.5));
            }
            else
            {
               str = Internat::ToDisplayString(lower);
            }
            item = new wxStaticText(mParent, 0, str);
            gridSizer->Add(item, 0, wxALIGN_CENTER_VERTICAL | wxALIGN_RIGHT | wxALL, 5);
         }
         else
         {
            gridSizer->Add(1, 1, 0);
         }

         mSliders[p] = new wxSlider(mParent, ID_SLIDERS + p,
                                    0, 0, 1000,
                                    wxDefaultPosition,
                                    wxSize(200, -1));
         mSliders[p]->SetName(labelText);
         gridSizer->Add(mSliders[p], 0, wxALIGN_CENTER_VERTICAL | wxEXPAND | wxALL, 5);
      
         if (hashi)
         {
            if (LADSPA_IS_HINT_INTEGER(hint.HintDescriptor) || forceint)
            {
               str.Printf(wxT("%d"), (int)(upper + 0.5));
            }
            else
            {
               str = Internat::ToDisplayString(upper);
            }
            item = new wxStaticText(mParent, 0, str);
            gridSizer->Add(item, 0, wxALIGN_CENTER_VERTICAL | wxALIGN_LEFT | wxALL, 5);
         }
         else
         {
            gridSizer->Add(1, 1, 0);
         }

         if (LADSPA_IS_HINT_INTEGER(hint.HintDescriptor) || forceint)
         {
            fieldText.Printf(wxT("%d"), (int)(mInputControls[p] + 0.5));

            wxIntegerValidator<float> vld(&mInputControls[p]);
            vld.SetRange(haslo ? lower : INT_MIN,
                           hashi ? upper : INT_MAX);
            mFields[p]->SetValidator(vld);
         }
         else
         {
            fieldText = Internat::ToDisplayString(mInputControls[p]);

            // > 12 decimal places can cause rounding errors in display.
            wxFloatingPointValidator<float> vld(12, &mInputControls[p]);
            vld.SetRange(haslo ? lower : -FLT_MAX,
                           hashi ? upper : FLT_MAX);

            // Set number of decimal places
            if (upper - lower < 10.0)
            {
               vld.SetStyle(wxNUM_VAL_THREE_TRAILING_ZEROES);
            }
            else if (upper - lower < 100.0)
            {
               vld.SetStyle(wxNUM_VAL_TWO_TRAILING_ZEROES);
            }
            else
            {
               vld.SetStyle(wxNUM_VAL_ONE_TRAILING_ZERO);
            }

            mFields[p]->SetValidator(vld);
         }

         // Set the textctrl value.  This will trigger an event so OnTextCtrl()
         // can update the slider.
         mFields[p]->SetValue(fieldText);
      }

      paramSizer->Add(gridSizer, 1, wxEXPAND | wxALL, 5);
      marginSizer->Add(paramSizer, 1, wxEXPAND | wxALL, 5);
   }

   if (mNumOutputControls > 0 )
   {
      wxSizer *paramSizer = new wxStaticBoxSizer(wxVERTICAL, mParent, _("Effect Output"));

      wxFlexGridSizer *gridSizer = new wxFlexGridSizer(2, 0, 0);
      gridSizer->AddGrowableCol(3);

      wxControl *item;

      for (unsigned long p = 0; p < mData->PortCount; p++)
      {
         LADSPA_PortDescriptor d = mData->PortDescriptors[p];
         if (LADSPA_IS_PORT_AUDIO(d) || LADSPA_IS_PORT_INPUT(d))
         {
            continue;
         }
         
         wxString labelText = LAT1CTOWX(mData->PortNames[p]);
         item = new wxStaticText(mParent, 0, labelText + wxT(":"));
         gridSizer->Add(item, 0, wxALIGN_CENTER_VERTICAL | wxALIGN_RIGHT | wxALL, 5);

         wxString fieldText;

         mFields[p] = new wxTextCtrl(mParent, wxID_ANY,
                                     fieldText,
                                     wxDefaultPosition,
                                     wxDefaultSize,
                                     wxTE_READONLY);
         mFields[p]->SetName(labelText);
         gridSizer->Add(mFields[p], 0, wxALIGN_CENTER_VERTICAL | wxALL, 5);
      }

      paramSizer->Add(gridSizer, 1, wxEXPAND | wxALL, 5);
      marginSizer->Add(paramSizer, 1, wxEXPAND | wxALL, 5);

      RefreshControls(true);
   }

   mParent->SetSizer(marginSizer);

   return true;
}

bool LadspaEffect::ValidateUI()
{
   if (!mParent->Validate())
   {
      return false;
   }

   if (GetType() == EffectTypeGenerate)
   {
      mHost->SetDuration(mDuration->GetValue());
   }

   return true;
}

bool LadspaEffect::HideUI()
{
   if (GetType() == EffectTypeAnalyze || mNumOutputControls > 0)
   {
      return false;
   }

   return true;
}

bool LadspaEffect::CloseUI()
{
   mParent->RemoveEventHandler(mEventHelper);
   delete mEventHelper;

   mUIHost = NULL;
   mParent = NULL;
   mDialog = NULL;

   return true;
}

void LadspaEffect::LoadUserPreset(const wxString & name)
{
   LoadParameters(name);
}

void LadspaEffect::SaveUserPreset(const wxString & name)
{
   SaveParameters(name);
}

void LadspaEffect::LoadFactoryPreset(int WXUNUSED(id))
{
   return;
}

void LadspaEffect::LoadFactoryDefaults()
{
   LoadParameters(mHost->GetFactoryDefaultsGroup());
}

wxArrayString LadspaEffect::GetFactoryPresets()
{
   return wxArrayString();
}

void LadspaEffect::ExportPresets()
{
}

void LadspaEffect::ImportPresets()
{
}

void LadspaEffect::ShowOptions()
{
}

// ============================================================================
// LadspaEffect Implementation
// ============================================================================

bool LadspaEffect::Load()
{
   if (mLib.IsLoaded())
   {
      return true;
   }

   LADSPA_Descriptor_Function mainFn = NULL;
   if (mLib.Load(mPath, wxDL_NOW))
   {
      wxLogNull logNo;

      mainFn = (LADSPA_Descriptor_Function) mLib.GetSymbol(wxT("ladspa_descriptor"));
      if (mainFn)
      {
         mData = mainFn(mIndex);
         return true;
      }
   }

   if (mLib.IsLoaded())
   {
      mLib.Unload();
   }

   return false;
}

void LadspaEffect::Unload()
{
   if (mLib.IsLoaded())
   {
      mLib.Unload();
   }
}

void LadspaEffect::LoadParameters(const wxString & group)
{
   wxString value;

   if (mHost->GetPrivateConfig(group, wxT("Value"), value, wxEmptyString))
   {
      wxStringTokenizer st(value, wxT(','));
      if (st.CountTokens() == mData->PortCount)
      {
         for (unsigned long p = 0; st.HasMoreTokens(); p++)
         {
            double val = 0.0;
            st.GetNextToken().ToDouble(&val);
            mInputControls[p] = val;
         }
      }
   }
}

void LadspaEffect::SaveParameters(const wxString & group)
{
   wxString parms;
   for (unsigned long p = 0; p < mData->PortCount; p++)
   {
      parms += wxString::Format(wxT(",%f"), mInputControls[p]);
   }

   mHost->SetPrivateConfig(group, wxT("Value"), parms.Mid(1));
}

LADSPA_Handle LadspaEffect::InitInstance(float sampleRate)
{
   /* Instantiate the plugin */
   LADSPA_Handle handle = mData->instantiate(mData, sampleRate);
   if (!handle)
   {
      return NULL;
   }

   for (unsigned long p = 0; p < mData->PortCount; p++)
   {
      LADSPA_PortDescriptor d = mData->PortDescriptors[p];
      if (LADSPA_IS_PORT_CONTROL(d))
      {
         if (LADSPA_IS_PORT_INPUT(d))
         {
            mData->connect_port(handle, p, &mInputControls[p]);
         }
         else
         {
            mData->connect_port(handle, p, &mOutputControls[p]);
         }
      }
   }

   if (mData->activate)
   {
      mData->activate(handle);
   }

   return handle;
}

void LadspaEffect::FreeInstance(LADSPA_Handle handle)
{
   if (mData->deactivate)
   {
      mData->deactivate(handle);
   }

   mData->cleanup(handle);
}

void LadspaEffect::OnCheckBox(wxCommandEvent & evt)
{
   int p = evt.GetId() - ID_TOGGLES;

   mInputControls[p] = mToggles[p]->GetValue();
}

void LadspaEffect::OnSlider(wxCommandEvent & evt)
{
   int p = evt.GetId() - ID_SLIDERS;

   float val;
   float lower = float(0.0);
   float upper = float(10.0);
   float range;
   bool forceint = false;

   LADSPA_PortRangeHint hint = mData->PortRangeHints[p];
   if (LADSPA_IS_HINT_BOUNDED_BELOW(hint.HintDescriptor))
      lower = hint.LowerBound;
   if (LADSPA_IS_HINT_BOUNDED_ABOVE(hint.HintDescriptor))
      upper = hint.UpperBound;
   if (LADSPA_IS_HINT_SAMPLE_RATE(hint.HintDescriptor)) {
      lower *= mSampleRate;
      upper *= mSampleRate;
      forceint = true;
   }

   range = upper - lower;

   val = (mSliders[p]->GetValue() / 1000.0) * range + lower;

   wxString str;
   if (LADSPA_IS_HINT_INTEGER(hint.HintDescriptor) || forceint)
      str.Printf(wxT("%d"), (int)(val + 0.5));
   else
      str = Internat::ToDisplayString(val);

   mFields[p]->SetValue(str);

   mInputControls[p] = val;
}

void LadspaEffect::OnTextCtrl(wxCommandEvent & evt)
{
   LadspaEffect *that = reinterpret_cast<LadspaEffect *>(this);
   int p = evt.GetId() - ID_TEXTS;

   float val;
   float lower = float(0.0);
   float upper = float(10.0);
   float range;

   val = Internat::CompatibleToDouble(that->mFields[p]->GetValue());

   LADSPA_PortRangeHint hint = that->mData->PortRangeHints[p];
   if (LADSPA_IS_HINT_BOUNDED_BELOW(hint.HintDescriptor))
      lower = hint.LowerBound;
   if (LADSPA_IS_HINT_BOUNDED_ABOVE(hint.HintDescriptor))
      upper = hint.UpperBound;
   if (LADSPA_IS_HINT_SAMPLE_RATE(hint.HintDescriptor)) {
      lower *= mSampleRate;
      upper *= mSampleRate;
   }
   range = upper - lower;

   if (val < lower)
      val = lower;
   if (val > upper)
      val = upper;

   mInputControls[p] = val;

   that->mSliders[p]->SetValue((int)(((val-lower)/range) * 1000.0 + 0.5));
}

void LadspaEffect::RefreshControls(bool outputOnly)
{
   if (!mParent)
   {
      return;
   }

   for (unsigned long p = 0; p < mData->PortCount; p++)
   {
      LADSPA_PortDescriptor d = mData->PortDescriptors[p];
      if (!(LADSPA_IS_PORT_CONTROL(d)))
      {
         continue;
      }

      wxString fieldText;
      LADSPA_PortRangeHint hint = mData->PortRangeHints[p];

      bool forceint = false;
      if (LADSPA_IS_HINT_SAMPLE_RATE(hint.HintDescriptor))
      {
         forceint = true;
      }

      if (LADSPA_IS_PORT_OUTPUT(d)) 
      {
         if (LADSPA_IS_HINT_INTEGER(hint.HintDescriptor) || forceint)
         {
            fieldText.Printf(wxT("%d"), (int)(mOutputControls[p] + 0.5));
         }
         else
         {
            fieldText = Internat::ToDisplayString(mOutputControls[p]);
         }

         mFields[p]->SetValue(fieldText);

         continue;
      }

      if (outputOnly)
      {
         continue;
      }

      if (LADSPA_IS_HINT_TOGGLED(hint.HintDescriptor))
      {
         mToggles[p]->SetValue(mInputControls[p] > 0);
         continue;
      }

      if (LADSPA_IS_HINT_INTEGER(hint.HintDescriptor) || forceint)
      {
         fieldText.Printf(wxT("%d"), (int)(mInputControls[p] + 0.5));
      }
      else
      {
         fieldText = Internat::ToDisplayString(mInputControls[p]);
      }

      // Set the textctrl value.  This will trigger an event so OnTextCtrl()
      // can update the slider.
      mFields[p]->SetValue(fieldText);
   }
}
