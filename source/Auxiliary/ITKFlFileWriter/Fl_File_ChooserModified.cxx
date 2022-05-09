// generated by Fast Light User Interface Designer (fluid) version 1.0104

#include "Fl_File_ChooserModified.H"

inline void Fl_File_ChooserModified::cb_window_i(Fl_Double_Window*, void*) {
  fileName->value("");
fileList->deselect();
Fl::remove_timeout((Fl_Timeout_Handler)previewCB, this);
window->hide();
}
void Fl_File_ChooserModified::cb_window(Fl_Double_Window* o, void* v) {
  ((Fl_File_ChooserModified*)(o->user_data()))->cb_window_i(o,v);
}

inline void Fl_File_ChooserModified::cb_showChoice_i(Fl_Choice*, void*) {
  showChoiceCB();
}
void Fl_File_ChooserModified::cb_showChoice(Fl_Choice* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->parent()->user_data()))->cb_showChoice_i(o,v);
}

inline void Fl_File_ChooserModified::cb_favoritesButton_i(Fl_Menu_Button*, void*) {
  favoritesButtonCB();
}
void Fl_File_ChooserModified::cb_favoritesButton(Fl_Menu_Button* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->parent()->user_data()))->cb_favoritesButton_i(o,v);
}

inline void Fl_File_ChooserModified::cb_newButton_i(Fl_Button*, void*) {
  newdir();
}
void Fl_File_ChooserModified::cb_newButton(Fl_Button* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->parent()->user_data()))->cb_newButton_i(o,v);
}

inline void Fl_File_ChooserModified::cb__i(Fl_Tile*, void*) {
  update_preview();
}
void Fl_File_ChooserModified::cb_(Fl_Tile* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->user_data()))->cb__i(o,v);
}

inline void Fl_File_ChooserModified::cb_fileList_i(Fl_File_Browser*, void*) {
  fileListCB();
}
void Fl_File_ChooserModified::cb_fileList(Fl_File_Browser* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->parent()->user_data()))->cb_fileList_i(o,v);
}

inline void Fl_File_ChooserModified::cb_previewButton_i(Fl_Check_Button*, void*) {
  preview(previewButton->value());
}
void Fl_File_ChooserModified::cb_previewButton(Fl_Check_Button* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->parent()->parent()->user_data()))->cb_previewButton_i(o,v);
}

inline void Fl_File_ChooserModified::cb_fileName_i(Fl_File_Input*, void*) {
  fileNameCB();
}
void Fl_File_ChooserModified::cb_fileName(Fl_File_Input* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->parent()->user_data()))->cb_fileName_i(o,v);
}

inline void Fl_File_ChooserModified::cb_okButton_i(Fl_Return_Button*, void*) {
  // Do any callback that is registered...
if (callback_)
  (*callback_)(this, data_);

window->hide();
}
void Fl_File_ChooserModified::cb_okButton(Fl_Return_Button* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->parent()->parent()->user_data()))->cb_okButton_i(o,v);
}

inline void Fl_File_ChooserModified::cb_Cancel_i(Fl_Button*, void*) {
  fileName->value("");
fileList->deselect();
Fl::remove_timeout((Fl_Timeout_Handler)previewCB, this);
window->hide();
}
void Fl_File_ChooserModified::cb_Cancel(Fl_Button* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->parent()->parent()->user_data()))->cb_Cancel_i(o,v);
}

inline void Fl_File_ChooserModified::cb_favList_i(Fl_File_Browser*, void*) {
  favoritesCB(favList);
}
void Fl_File_ChooserModified::cb_favList(Fl_File_Browser* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->user_data()))->cb_favList_i(o,v);
}

inline void Fl_File_ChooserModified::cb_favUpButton_i(Fl_Button*, void*) {
  favoritesCB(favUpButton);
}
void Fl_File_ChooserModified::cb_favUpButton(Fl_Button* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->user_data()))->cb_favUpButton_i(o,v);
}

inline void Fl_File_ChooserModified::cb_favDeleteButton_i(Fl_Button*, void*) {
  favoritesCB(favDeleteButton);
}
void Fl_File_ChooserModified::cb_favDeleteButton(Fl_Button* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->user_data()))->cb_favDeleteButton_i(o,v);
}

inline void Fl_File_ChooserModified::cb_favDownButton_i(Fl_Button*, void*) {
  favoritesCB(favDownButton);
}
void Fl_File_ChooserModified::cb_favDownButton(Fl_Button* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->user_data()))->cb_favDownButton_i(o,v);
}

inline void Fl_File_ChooserModified::cb_favCancelButton_i(Fl_Button*, void*) {
  favWindow->hide();
}
void Fl_File_ChooserModified::cb_favCancelButton(Fl_Button* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->user_data()))->cb_favCancelButton_i(o,v);
}

inline void Fl_File_ChooserModified::cb_favOkButton_i(Fl_Return_Button*, void*) {
  favoritesCB(favOkButton);
}
void Fl_File_ChooserModified::cb_favOkButton(Fl_Return_Button* o, void* v) {
  ((Fl_File_ChooserModified*)(o->parent()->user_data()))->cb_favOkButton_i(o,v);
}

Fl_File_ChooserModified::Fl_File_ChooserModified(const char *d, const char *p, int t, const char *title) {
  Fl_Double_Window* w;
  { Fl_Double_Window* o = window = new Fl_Double_Window(490, 380, "Choose File");
    w = o;
    o->callback((Fl_Callback*)cb_window, (void*)(this));
    { Fl_Group* o = new Fl_Group(65, 10, 415, 25);
      { Fl_Choice* o = showChoice = new Fl_Choice(65, 10, 215, 25, "Show:");
        o->down_box(FL_BORDER_BOX);
        o->callback((Fl_Callback*)cb_showChoice);
        Fl_Group::current()->resizable(o);
        showChoice->label(show_label);
      }
      { Fl_Menu_Button* o = favoritesButton = new Fl_Menu_Button(290, 10, 155, 25, "Favorites");
        o->down_box(FL_BORDER_BOX);
        o->callback((Fl_Callback*)cb_favoritesButton);
        o->align(FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
        favoritesButton->label(favorites_label);
      }
      { Fl_Button* o = newButton = new Fl_Button(455, 10, 25, 25);
        o->tooltip("Create a new directory.");
        o->labelsize(8);
        o->callback((Fl_Callback*)cb_newButton);
      }
      o->end();
    }
    { Fl_Tile* o = new Fl_Tile(10, 45, 470, 225);
      o->callback((Fl_Callback*)cb_);
      { Fl_File_Browser* o = fileList = new Fl_File_Browser(10, 45, 295, 225);
        o->type(2);
        o->callback((Fl_Callback*)cb_fileList);
        w->hotspot(o);
      }
      { Fl_Box* o = previewBox = new Fl_Box(305, 45, 175, 225, "?");
        o->box(FL_DOWN_BOX);
        o->labelsize(100);
        o->align(FL_ALIGN_CLIP|FL_ALIGN_INSIDE);
      }
      o->end();
      Fl_Group::current()->resizable(o);
    }
    { Fl_Group* o = new Fl_Group(0, 275, 480, 95);
      { Fl_Group* o = new Fl_Group(10, 275, 470, 20);
        { Fl_Check_Button* o = previewButton = new Fl_Check_Button(10, 275, 170, 20, "Preview");
          o->down_box(FL_DOWN_BOX);
          o->value(1);
          o->shortcut(0x80070);
          o->callback((Fl_Callback*)cb_previewButton);
          previewButton->label(preview_label);
        }
        { Fl_Box* o = new Fl_Box(10, 275, 395, 20);
          Fl_Group::current()->resizable(o);
        }
        o->end();
      }
      { Fl_File_Input* o = fileName = new Fl_File_Input(115, 300, 365, 35);
        o->callback((Fl_Callback*)cb_fileName);
        o->when(FL_WHEN_ENTER_KEY);
        Fl_Group::current()->resizable(o);
        fileName->when(FL_WHEN_CHANGED | FL_WHEN_ENTER_KEY_ALWAYS);
      }
      { Fl_Box* o = new Fl_Box(10, 310, 105, 25, "Filename:");
        o->align(FL_ALIGN_RIGHT|FL_ALIGN_INSIDE);
        o->label(filename_label);
      }
      { Fl_Group* o = new Fl_Group(10, 345, 470, 25);
        { Fl_Return_Button* o = okButton = new Fl_Return_Button(300, 345, 85, 25, "OK");
          o->callback((Fl_Callback*)cb_okButton);
          okButton->label(fl_ok);
        }
        { Fl_Button* o = new Fl_Button(395, 345, 85, 25, "Cancel");
          o->callback((Fl_Callback*)cb_Cancel);
          o->label(fl_cancel);
        }
        { Fl_Box* o = new Fl_Box(10, 345, 300, 25);
          Fl_Group::current()->resizable(o);
        }
        o->end();
      }
      o->end();
    }
    if (title) window->label(title);
    o->set_modal();
    o->end();
  }
  { Fl_Double_Window* o = favWindow = new Fl_Double_Window(355, 150, "Manage Favorites");
    w = o;
    o->user_data((void*)(this));
    { Fl_File_Browser* o = favList = new Fl_File_Browser(10, 10, 300, 95);
      o->type(2);
      o->callback((Fl_Callback*)cb_favList);
    }
    { Fl_Button* o = favUpButton = new Fl_Button(320, 10, 25, 25, "@8>");
      o->callback((Fl_Callback*)cb_favUpButton);
    }
    { Fl_Button* o = favDeleteButton = new Fl_Button(320, 45, 25, 25, "X");
      o->labelfont(1);
      o->callback((Fl_Callback*)cb_favDeleteButton);
    }
    { Fl_Button* o = favDownButton = new Fl_Button(320, 80, 25, 25, "@2>");
      o->callback((Fl_Callback*)cb_favDownButton);
    }
    { Fl_Button* o = favCancelButton = new Fl_Button(270, 115, 75, 25, "Cancel");
      o->callback((Fl_Callback*)cb_favCancelButton);
      favCancelButton->label(fl_cancel);
    }
    { Fl_Return_Button* o = favOkButton = new Fl_Return_Button(185, 115, 75, 25, "OK");
      o->callback((Fl_Callback*)cb_favOkButton);
      favOkButton->label(fl_ok);
    }
    favWindow->label(manage_favorites_label);
    o->set_modal();
    o->end();
  }
  callback_ = 0;
data_ = 0;
directory_[0] = 0;
window->size_range(window->w(), window->h(), Fl::w(), Fl::h());
type(t);
filter(p);
update_favorites();
value(d);
type(t);
int e;
prefs_.get("preview", e, 1);
preview(e);
}

Fl_File_ChooserModified::~Fl_File_ChooserModified() {
  Fl::remove_timeout((Fl_Timeout_Handler)previewCB, this);
delete window;
delete favWindow;
}

void Fl_File_ChooserModified::callback(void (*cb)(Fl_File_ChooserModified *, void *), void *d ) {
  callback_ = cb;
data_     = d;
}

void Fl_File_ChooserModified::color(Fl_Color c) {
  fileList->color(c);
}

Fl_Color Fl_File_ChooserModified::color() {
  return (fileList->color());
}

char * Fl_File_ChooserModified::directory() {
  return directory_;
}

const char * Fl_File_ChooserModified::filter() {
  return (fileList->filter());
}

int Fl_File_ChooserModified::filter_value() {
  return showChoice->value();
}

void Fl_File_ChooserModified::filter_value(int f) {
  showChoice->value(f);
showChoiceCB();
}

void Fl_File_ChooserModified::hide() {
  window->hide();
}

void Fl_File_ChooserModified::iconsize(uchar s) {
  fileList->iconsize(s);
}

uchar Fl_File_ChooserModified::iconsize() {
  return (fileList->iconsize());
}

void Fl_File_ChooserModified::label(const char *l) {
  window->label(l);
}

const char * Fl_File_ChooserModified::label() {
  return (window->label());
}

void Fl_File_ChooserModified::show() {
  window->hotspot(fileList);
window->show();
fileName->take_focus();
}

int Fl_File_ChooserModified::shown() {
  return window->shown();
}

void Fl_File_ChooserModified::textcolor(Fl_Color c) {
  fileList->textcolor(c);
}

Fl_Color Fl_File_ChooserModified::textcolor() {
  return (fileList->textcolor());
}

void Fl_File_ChooserModified::textfont(uchar f) {
  fileList->textfont(f);
}

uchar Fl_File_ChooserModified::textfont() {
  return (fileList->textfont());
}

void Fl_File_ChooserModified::textsize(uchar s) {
  fileList->textsize(s);
}

uchar Fl_File_ChooserModified::textsize() {
  return (fileList->textsize());
}

void Fl_File_ChooserModified::type(int t) {
  type_ = t;
if (t & MULTI)
  fileList->type(FL_MULTI_BROWSER);
else
  fileList->type(FL_HOLD_BROWSER);
if (t & CREATE)
  newButton->activate();
else
  newButton->deactivate();
if (t & DIRECTORY)
  fileList->filetype(Fl_File_Browser::DIRECTORIES);
else
  fileList->filetype(Fl_File_Browser::FILES);
}

int Fl_File_ChooserModified::type() {
  return (type_);
}

int Fl_File_ChooserModified::visible() {
  return window->visible();
}
