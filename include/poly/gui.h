#pragma once

#include <poly/common.h>
#include <poly/subscreen.h>
#include <nanogui/screen.h>

class GUI : public nanogui::Screen {
    GUI();

    virtual bool keyboardEvent(int key, int scancode, int action, int modifiers);
    virtual void draw(NVGcontext* ctx);
    virtual void drawContents();

    void initializeGUI();
};
