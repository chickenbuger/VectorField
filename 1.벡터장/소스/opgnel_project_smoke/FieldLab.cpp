#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <shellapi.h>
#include <wincodec.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include "GL/glut.h"
#include "FieldMath.h"
#include "../tests/FieldMathTests.h"

#pragma comment(lib, "windowscodecs.lib")
#pragma comment(lib, "ole32.lib")
#pragma comment(lib, "shell32.lib")

void LegacyMain(int argc, char** argv);

namespace fieldlab {
namespace {
constexpr double CanvasWidth = 1180.0;
constexpr double CanvasHeight = 860.0;
constexpr double PlotX = 54.0;
constexpr double PlotY = 128.0;
constexpr double PlotSize = 660.0;
constexpr double PanelX = 770.0;
constexpr double Domain = 2.0;

FieldType fieldType = FieldType::Rotation;
double strength = 1.0;
int gridSize = 21;
bool normalizedArrows = false;
Vec2 selected{1.0, 0.5};
int windowWidth = static_cast<int>(CanvasWidth);
int windowHeight = static_cast<int>(CanvasHeight);
double viewportScale = 1.0;
int viewportX = 0;
int viewportYTop = 0;
bool captureRequested = false;
std::wstring captureDirectory = L"captures";
std::string captureStatus = "F12 saves the rendered frame as PNG";

struct Color { double r, g, b; };
const Color Background{0.043, 0.063, 0.09};
const Color Surface{0.071, 0.098, 0.137};
const Color Muted{0.48, 0.58, 0.69};
const Color TextColor{0.87, 0.92, 0.97};
const Color Accent{0.32, 0.87, 0.77};

void SetColor(Color c) { glColor3d(c.r, c.g, c.b); }

void Rect(double x, double y, double w, double h, Color c) {
    SetColor(c);
    glBegin(GL_QUADS);
    glVertex2d(x, y); glVertex2d(x + w, y);
    glVertex2d(x + w, y + h); glVertex2d(x, y + h);
    glEnd();
}

void Line(double x1, double y1, double x2, double y2, Color c, float width = 1.0f) {
    SetColor(c);
    glLineWidth(width);
    glBegin(GL_LINES);
    glVertex2d(x1, y1); glVertex2d(x2, y2);
    glEnd();
    glLineWidth(1.0f);
}

void Text(double x, double y, const std::string& text, Color c = TextColor,
          void* font = GLUT_BITMAP_HELVETICA_12) {
    SetColor(c);
    glRasterPos2d(x, y);
    for (unsigned char ch : text) glutBitmapCharacter(font, ch);
}

std::string Number(double value, int digits = 3) {
    if (std::abs(value) < 0.5 * std::pow(10.0, -digits)) value = 0.0;
    char result[80];
    std::snprintf(result, sizeof(result), "%.*f", digits, value);
    return result;
}

const char* FieldName(FieldType type) {
    switch (type) {
    case FieldType::Uniform: return "Uniform field";
    case FieldType::Rotation: return "Rotational field";
    case FieldType::Source: return "Divergent field";
    case FieldType::Saddle: return "Saddle field";
    }
    return "";
}

const char* Formula(FieldType type) {
    switch (type) {
    case FieldType::Uniform: return "F(x,y) = k (1, 0.25)";
    case FieldType::Rotation: return "F(x,y) = k (-y, x)";
    case FieldType::Source: return "F(x,y) = k (x, y)";
    case FieldType::Saddle: return "F(x,y) = k (x, -y)";
    }
    return "";
}

Vec2 PlotPosition(Vec2 p) {
    return {PlotX + (p.x + Domain) / (2.0 * Domain) * PlotSize,
            PlotY + (Domain - p.y) / (2.0 * Domain) * PlotSize};
}

Color MagnitudeColor(double magnitude) {
    const double t = magnitude / (1.0 + magnitude);
    return {0.24 + 0.73 * t, 0.80 - 0.30 * t, 0.94 - 0.65 * t};
}

void DrawArrow(Vec2 position, Vec2 vector, Color color, double spacing,
               bool highlight = false) {
    const Vec2 base = PlotPosition(position);
    const double magnitude = vector.Length();
    if (magnitude < 1e-12) {
        Rect(base.x - 2, base.y - 2, 4, 4, color);
        return;
    }
    const Vec2 direction = Normalize(vector);
    const double fraction = normalizedArrows ? 0.72 : 0.80 * magnitude / (1.0 + magnitude);
    const double length = spacing * fraction;
    const Vec2 delta{direction.x * length, -direction.y * length};
    const Vec2 tip = base + delta;
    const Vec2 side{-delta.y / length, delta.x / length};
    const double head = std::min(highlight ? 8.0 : 5.0, length * 0.38);
    const Vec2 back = tip - delta * (head / length);
    Line(base.x, base.y, tip.x, tip.y, color, highlight ? 2.5f : 1.3f);
    Line(tip.x, tip.y, back.x + side.x * head * 0.50,
         back.y + side.y * head * 0.50, color, highlight ? 2.0f : 1.3f);
    Line(tip.x, tip.y, back.x - side.x * head * 0.50,
         back.y - side.y * head * 0.50, color, highlight ? 2.0f : 1.3f);
}

void DrawPlot() {
    Rect(PlotX, PlotY, PlotSize, PlotSize, {0.055, 0.078, 0.109});
    for (int i = 0; i <= 8; ++i) {
        const double offset = PlotSize * i / 8.0;
        const Color c = i == 4 ? Color{0.29, 0.36, 0.44} : Color{0.10, 0.14, 0.19};
        Line(PlotX + offset, PlotY, PlotX + offset, PlotY + PlotSize, c);
        Line(PlotX, PlotY + offset, PlotX + PlotSize, PlotY + offset, c);
        if (i % 2 == 0) {
            Text(PlotX + offset - 7, PlotY + PlotSize + 22, Number(-2 + i * 0.5, 0), Muted);
            Text(PlotX - 26, PlotY + offset + 4, Number(2 - i * 0.5, 0), Muted);
        }
    }
    Text(PlotX + PlotSize - 12, PlotY + PlotSize + 43, "x", Muted);
    Text(PlotX - 25, PlotY - 14, "y", Muted);

    const double step = 2.0 * Domain / (gridSize - 1);
    const double pixelSpacing = PlotSize / (gridSize - 1);
    // Samples include the domain boundary; clip arrows to the plot rectangle.
    glEnable(GL_SCISSOR_TEST);
    glScissor(viewportX + static_cast<int>(PlotX * viewportScale),
              windowHeight - viewportYTop - static_cast<int>((PlotY + PlotSize) * viewportScale),
              static_cast<int>(PlotSize * viewportScale),
              static_cast<int>(PlotSize * viewportScale));
    for (int j = 0; j < gridSize; ++j) {
        for (int i = 0; i < gridSize; ++i) {
            Vec2 p{-Domain + step * i, -Domain + step * j};
            Vec2 vector = Evaluate(fieldType, p, strength);
            DrawArrow(p, vector, MagnitudeColor(vector.Length()), pixelSpacing);
        }
    }

    const Vec2 pixel = PlotPosition(selected);
    Line(pixel.x - 9, pixel.y, pixel.x + 9, pixel.y, TextColor, 1.5f);
    Line(pixel.x, pixel.y - 9, pixel.x, pixel.y + 9, TextColor, 1.5f);
    DrawArrow(selected, Evaluate(fieldType, selected, strength), TextColor, 72.0, true);
    glDisable(GL_SCISSOR_TEST);
}

void DrawPanel() {
    Rect(PanelX, 110, 364, 690, Surface);
    Text(PanelX + 22, 143, "FIELD PRESETS", Muted);
    const FieldType types[] = {FieldType::Uniform, FieldType::Rotation,
                               FieldType::Source, FieldType::Saddle};
    for (int i = 0; i < 4; ++i) {
        const double y = 158 + i * 47.0;
        const bool active = types[i] == fieldType;
        Rect(PanelX + 16, y, 332, 38, active ? Color{0.11, 0.24, 0.25} : Color{0.086, 0.12, 0.16});
        if (active) Rect(PanelX + 16, y, 3, 38, Accent);
        Text(PanelX + 29, y + 24, std::to_string(i + 1) + "   " + FieldName(types[i]),
             active ? Accent : TextColor, GLUT_BITMAP_HELVETICA_18);
    }
    Text(PanelX + 22, 384, Formula(fieldType), Accent, GLUT_BITMAP_HELVETICA_18);
    Text(PanelX + 22, 414, "Strength k                 " + Number(strength, 2));
    Text(PanelX + 22, 440, "Grid samples              " + std::to_string(gridSize) + " x " + std::to_string(gridSize));
    Text(PanelX + 22, 466, normalizedArrows ? "Arrow length               normalized" : "Arrow length               compressed magnitude");
    Line(PanelX + 22, 488, PanelX + 342, 488, {0.16, 0.21, 0.27});
    Text(PanelX + 22, 516, "PROBE / CLICK IN THE FIELD", Muted);
    const Vec2 value = Evaluate(fieldType, selected, strength);
    const Differential differential = Derivatives(fieldType, selected, strength);
    Text(PanelX + 22, 545, "p = (" + Number(selected.x) + ", " + Number(selected.y) + ")");
    Text(PanelX + 22, 571, "F = (" + Number(value.x) + ", " + Number(value.y) + ")", Accent);
    Text(PanelX + 22, 597, "|F| = " + Number(value.Length()));
    Text(PanelX + 22, 633, "Divergence       " + Number(differential.divergence));
    Text(PanelX + 22, 659, "Curl (z)             " + Number(differential.curl));
    Text(PanelX + 22, 687, "Central differences, h = 0.0001", Muted);
    Line(PanelX + 22, 708, PanelX + 342, 708, {0.16, 0.21, 0.27});
    Text(PanelX + 22, 738, "1-4  preset     +/-  strength     G  grid", Muted);
    Text(PanelX + 22, 762, "N  normalize     R  reset     F12  capture", Muted);
    Text(PanelX + 22, 786, "Click the preset buttons or the plot", Muted);
}

std::wstring CaptureName() {
    switch (fieldType) {
    case FieldType::Rotation: return L"01_\ud68c\uc804\uc7a5.png";
    case FieldType::Source: return L"02_\ubc1c\uc0b0\uc7a5.png";
    case FieldType::Saddle: return L"03_\uc548\uc7a5\uc7a5.png";
    case FieldType::Uniform: return L"04_\uade0\uc77c\uc7a5.png";
    }
    return L"capture.png";
}

// Save pixels from the real OpenGL back buffer, including the panel and probe.
// WIC writes PNG without adding another third-party graphics dependency.
bool SaveFrame() {
    if (!CreateDirectoryW(captureDirectory.c_str(), nullptr) && GetLastError() != ERROR_ALREADY_EXISTS)
        return false;
    const std::wstring path = captureDirectory + L"\\" + CaptureName();
    const UINT width = static_cast<UINT>(windowWidth);
    const UINT height = static_cast<UINT>(windowHeight);
    std::vector<BYTE> rgb(static_cast<size_t>(width) * height * 3);
    std::vector<BYTE> bgr(rgb.size());
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadBuffer(GL_BACK);
    glFinish();
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, rgb.data());
    for (UINT y = 0; y < height; ++y) {
        for (UINT x = 0; x < width; ++x) {
            const size_t out = (static_cast<size_t>(y) * width + x) * 3;
            const size_t in = (static_cast<size_t>(height - 1 - y) * width + x) * 3;
            bgr[out] = rgb[in + 2]; bgr[out + 1] = rgb[in + 1]; bgr[out + 2] = rgb[in];
        }
    }
    const HRESULT init = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);
    if (FAILED(init) && init != RPC_E_CHANGED_MODE) return false;
    IWICImagingFactory* factory = nullptr;
    IWICStream* stream = nullptr;
    IWICBitmapEncoder* encoder = nullptr;
    IWICBitmapFrameEncode* frame = nullptr;
    HRESULT hr = CoCreateInstance(CLSID_WICImagingFactory, nullptr, CLSCTX_INPROC_SERVER,
                                  IID_PPV_ARGS(&factory));
    if (SUCCEEDED(hr)) hr = factory->CreateStream(&stream);
    if (SUCCEEDED(hr)) hr = stream->InitializeFromFilename(path.c_str(), GENERIC_WRITE);
    if (SUCCEEDED(hr)) hr = factory->CreateEncoder(GUID_ContainerFormatPng, nullptr, &encoder);
    if (SUCCEEDED(hr)) hr = encoder->Initialize(stream, WICBitmapEncoderNoCache);
    if (SUCCEEDED(hr)) hr = encoder->CreateNewFrame(&frame, nullptr);
    if (SUCCEEDED(hr)) hr = frame->Initialize(nullptr);
    if (SUCCEEDED(hr)) hr = frame->SetSize(width, height);
    WICPixelFormatGUID format = GUID_WICPixelFormat24bppBGR;
    if (SUCCEEDED(hr)) hr = frame->SetPixelFormat(&format);
    if (SUCCEEDED(hr) && !IsEqualGUID(format, GUID_WICPixelFormat24bppBGR)) hr = E_FAIL;
    if (SUCCEEDED(hr)) hr = frame->WritePixels(height, width * 3, static_cast<UINT>(bgr.size()), bgr.data());
    if (SUCCEEDED(hr)) hr = frame->Commit();
    if (SUCCEEDED(hr)) hr = encoder->Commit();
    if (frame) frame->Release();
    if (encoder) encoder->Release();
    if (stream) stream->Release();
    if (factory) factory->Release();
    if (SUCCEEDED(init)) CoUninitialize();
    return SUCCEEDED(hr);
}

void Display() {
    glClearColor(static_cast<GLfloat>(Background.r), static_cast<GLfloat>(Background.g),
                 static_cast<GLfloat>(Background.b), 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    Text(54, 42, "VECTORFIELD  /  DEVELOPMENT LAB", Accent);
    Text(54, 77, "01    Vector fields", TextColor, GLUT_BITMAP_HELVETICA_18);
    Text(54, 102, "Position -> vector -> magnitude -> divergence & curl", Muted);
    Text(PanelX + 165, 45, "STAGE 01 / 06", Muted);
    DrawPlot();
    DrawPanel();
    Text(54, 842, "Color: low magnitude (cyan) -> high (amber).  Display length does not change F.", Muted);
    Text(PanelX, 829, captureStatus, Muted);
    if (captureRequested) {
        captureRequested = false;
        const bool saved = SaveFrame();
        captureStatus = saved ? "PNG saved. F12 replaces this preset capture." : "Capture failed. Check the output folder.";
        std::printf("Capture %s (preset %s)\n", saved ? "saved" : "FAILED", FieldName(fieldType));
        std::fflush(stdout);
        glutPostRedisplay();
    }
    glutSwapBuffers();
}

void Reshape(int width, int height) {
    windowWidth = std::max(width, 1);
    windowHeight = std::max(height, 1);
    viewportScale = std::min(windowWidth / CanvasWidth, windowHeight / CanvasHeight);
    const int viewWidth = static_cast<int>(CanvasWidth * viewportScale);
    const int viewHeight = static_cast<int>(CanvasHeight * viewportScale);
    viewportX = (windowWidth - viewWidth) / 2;
    viewportYTop = (windowHeight - viewHeight) / 2;
    glViewport(viewportX, windowHeight - viewportYTop - viewHeight, viewWidth, viewHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, CanvasWidth, CanvasHeight, 0.0, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glutPostRedisplay();
}

void Keyboard(unsigned char key, int, int) {
    switch (key) {
    case '1': fieldType = FieldType::Uniform; break;
    case '2': fieldType = FieldType::Rotation; break;
    case '3': fieldType = FieldType::Source; break;
    case '4': fieldType = FieldType::Saddle; break;
    case '+': case '=': strength = std::min(3.0, strength + 0.25); break;
    case '-': case '_': strength = std::max(0.25, strength - 0.25); break;
    case 'g': case 'G': gridSize = gridSize == 21 ? 31 : (gridSize == 31 ? 41 : 21); break;
    case 'n': case 'N': normalizedArrows = !normalizedArrows; break;
    case 'r': case 'R':
        fieldType = FieldType::Rotation; strength = 1.0; gridSize = 21;
        normalizedArrows = false; selected = {1.0, 0.5}; break;
    }
    glutPostRedisplay();
}

void Special(int key, int, int) {
    if (key == GLUT_KEY_F12) { captureRequested = true; glutPostRedisplay(); }
}

void Mouse(int button, int state, int x, int y) {
    if (button != GLUT_LEFT_BUTTON || state != GLUT_DOWN) return;
    const double px = (x - viewportX) / viewportScale;
    const double py = (y - viewportYTop) / viewportScale;
    if (px >= PlotX && px <= PlotX + PlotSize && py >= PlotY && py <= PlotY + PlotSize) {
        selected = {(px - PlotX) / PlotSize * 4.0 - 2.0,
                     2.0 - (py - PlotY) / PlotSize * 4.0};
    } else if (px >= PanelX + 16 && px <= PanelX + 348) {
        for (int i = 0; i < 4; ++i) {
            if (py >= 158 + i * 47.0 && py <= 196 + i * 47.0)
                fieldType = static_cast<FieldType>(i);
        }
    }
    glutPostRedisplay();
}
} // namespace

int Run(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--self-test") == 0) return RunFieldMathTests();
        if (std::strcmp(argv[i], "--legacy") == 0) { LegacyMain(argc, argv); return 0; }
    }
    int wideCount = 0;
    LPWSTR* wideArgs = CommandLineToArgvW(GetCommandLineW(), &wideCount);
    if (wideArgs) {
        for (int i = 1; i + 1 < wideCount; ++i) {
            if (std::wstring(wideArgs[i]) == L"--capture-dir") captureDirectory = wideArgs[++i];
        }
        LocalFree(wideArgs);
    }
    // Keep custom Unicode arguments away from GLUT's legacy argument parser.
    int glutCount = 1;
    glutInit(&glutCount, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(30, 35);
    glutCreateWindow("VectorField - 01 Vector fields");
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    glutKeyboardFunc(Keyboard);
    glutSpecialFunc(Special);
    glutMouseFunc(Mouse);
    // Event-driven drawing avoids an idle loop consuming a CPU core.
    glutMainLoop();
    return 0;
}
} // namespace fieldlab

int main(int argc, char** argv) { return fieldlab::Run(argc, argv); }
