#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <shellapi.h>
#include <wincodec.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <string>
#include <vector>
#include "GL/glut.h"
#include "FieldMath.h"
#include "DynamicsMath.h"
#include "FluidGrid.h"
#include "../tests/FieldMathTests.h"
#include "../tests/DynamicsTests.h"
#include "../tests/FluidGridTests.h"

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
int stage = 1;
constexpr int AvailableStages = 5;
FluidGrid fluid(64);
bool projectionEnabled = true;
bool buoyancyEnabled = true;
bool emitterEnabled = true;
bool mouseHeld = false;
Vec2 lastMousePosition;
double viscosity = 0.00015;
double diffusion = 0.00002;
bool requestedNoProjection = false;
struct Tracer { Vec2 position; std::deque<Vec2> trail; };
std::vector<Tracer> tracers;
bool paused = false;
double accumulator = 0.0;
std::vector<GravitySource> gravitySources;
int selectedSource = 0;
bool draggingSource = false;
double simulationTime = 0.0;
double captureAfter = -1.0;
bool quitAfterCapture = false;
std::wstring captureNameOverride;
int lastTimerTime = 0;
double recordElapsed = 0.0;
double requestedStrength = -1.0;
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

void ResetStage() {
    fieldType = FieldType::Rotation; strength = 1.0; gridSize = 21;
    normalizedArrows = false; selected = {1.0, 0.5}; simulationTime = 0.0;
    gravitySources = {{{-0.65, 0.2}, 1.0}, {{0.8, -0.1}, 1.5}};
    selectedSource = 0; draggingSource = false;
    mouseHeld = false; projectionEnabled = true; buoyancyEnabled = true; emitterEnabled = true;
    viscosity = 0.00015; diffusion = stage == 4 ? 0.001 : 0.00002;
    fluid.SeedDemo();
    if (stage == 4) strength = 0.65;
    if (stage == 5) fluid.AddVelocity(0.52,0.5,1.2,0.5,0.16);
    if (stage == 6) fluid.Clear();
    paused = false; accumulator = 0.0; recordElapsed = 0.0;
    tracers.clear();
    for (int j = 0; j < 18; ++j) for (int i = 0; i < 18; ++i)
        tracers.push_back({{-1.4 + 2.8 * i / 17.0, -1.4 + 2.8 * j / 17.0}, {}});
    const char* titles[] = {"", "VectorField - 01 Vector fields", "VectorField - 02 Gravity field", "VectorField - 03 Particle advection", "VectorField - 04 Advection and diffusion", "VectorField - 05 Incompressible fluid", "VectorField - 06 Interactive smoke"};
    glutSetWindowTitle(titles[stage]);
}

Vec2 CurrentField(Vec2 position) {
    if (stage >= 4) return fluid.SampleVelocity((position + Vec2{2,2}) / 4.0);
    return stage == 2 ? GravityField(position, gravitySources, 0.4 * strength, 0.12)
                      : Evaluate(fieldType, position, strength);
}

Differential CurrentDerivatives(Vec2 position) {
    if (stage != 2) return Derivatives(fieldType, position, strength);
    const double h = 1e-4;
    const Vec2 dx = (CurrentField(position + Vec2{h, 0}) - CurrentField(position - Vec2{h, 0})) / (2 * h);
    const Vec2 dy = (CurrentField(position + Vec2{0, h}) - CurrentField(position - Vec2{0, h})) / (2 * h);
    return {dx.x + dy.y, dx.y - dy.x};
}

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

std::string Scientific(double value) {
    char result[48]; std::snprintf(result,sizeof(result),"%.3e",value); return result;
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

    if (stage >= 4) {
        const double cell = PlotSize / fluid.Size();
        for (int j = 1; j <= fluid.Size(); ++j) for (int i = 1; i <= fluid.Size(); ++i) {
            const double density = fluid.Dye(i,j);
            const double intensity = 1.0 - std::exp(-1.3 * density);
            const Color color = stage == 6 ? Color{0.10+0.76*intensity,0.12+0.78*intensity,0.16+0.82*intensity}
                : Color{0.06+0.11*intensity,0.09+0.69*intensity,0.12+0.85*intensity};
            Rect(PlotX+(i-1)*cell, PlotY+(fluid.Size()-j)*cell,cell+0.05,cell+0.05,color);
        }
    }

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
            Vec2 vector = CurrentField(p);
            DrawArrow(p, vector, stage >= 4 ? Color{0.50,0.59,0.48} : MagnitudeColor(vector.Length()), pixelSpacing);
        }
    }

    const Vec2 pixel = PlotPosition(selected);
    Line(pixel.x - 9, pixel.y, pixel.x + 9, pixel.y, TextColor, 1.5f);
    Line(pixel.x, pixel.y - 9, pixel.x, pixel.y + 9, TextColor, 1.5f);
    DrawArrow(selected, CurrentField(selected), TextColor, 72.0, true);
    if (stage == 2) {
        for (size_t i = 0; i < gravitySources.size(); ++i) {
            const Vec2 center = PlotPosition(gravitySources[i].position);
            const double radius = 11.0 + gravitySources[i].mass * 3.0;
            SetColor(i == static_cast<size_t>(selectedSource) ? Accent : Color{0.98, 0.70, 0.28});
            glBegin(GL_TRIANGLE_FAN);
            glVertex2d(center.x, center.y);
            for (int a = 0; a <= 32; ++a) {
                const double angle = a * 6.283185307179586 / 32.0;
                glVertex2d(center.x + radius * std::cos(angle), center.y + radius * std::sin(angle));
            }
            glEnd();
            Text(center.x + radius + 4, center.y - 6, "m=" + Number(gravitySources[i].mass, 1));
        }
    }
    if (stage == 3) {
        for (const auto& tracer : tracers) {
            glColor4d(0.30, 0.82, 0.94, 0.55);
            glBegin(GL_LINE_STRIP);
            for (const auto& p : tracer.trail) { const Vec2 pixel = PlotPosition(p); glVertex2d(pixel.x,pixel.y); }
            glEnd();
        }
        SetColor(Accent); glPointSize(4.0f); glBegin(GL_POINTS);
        for (const auto& tracer : tracers) { const Vec2 pixel = PlotPosition(tracer.position); glVertex2d(pixel.x,pixel.y); }
        glEnd(); glPointSize(1.0f);
    }
    glDisable(GL_SCISSOR_TEST);
}

void DrawFluidPanel() {
    Rect(PanelX,110,364,690,Surface);
    Text(PanelX+22,143,"EULERIAN / MAC GRID SOLVER",Muted);
    const char* mode = stage == 4 ? "Passive dye transport" : stage == 5 ? "Pressure projection" : "Interactive smoke";
    Text(PanelX+22,184,mode,Accent,GLUT_BITMAP_HELVETICA_18);
    Text(PanelX+22,218,stage == 4 ? "Given velocity + dye advection" : "Velocity advection + viscosity");
    Text(PanelX+22,251,stage == 4 ? "Projection: not part of this stage" : (projectionEnabled ? "Projection: ON / SOR 100" : "Projection: OFF"),stage == 4 ? Muted : Accent);
    Text(PanelX+22,284,stage == 6 ? (buoyancyEnabled ? "Buoyancy: ON" : "Buoyancy: OFF") : "Buoyancy: OFF",Muted);
    Text(PanelX+22,317,stage == 6 ? (emitterEnabled ? "Two emitters: ON" : "Two emitters: OFF") : stage == 4 ? std::string(Formula(fieldType)) + "  k=" + Number(strength,2) : "Velocity: simulated, not prescribed",Muted);
    Line(PanelX+22,339,PanelX+342,339,{0.16,0.21,0.27});
    Text(PanelX+22,369,"t = " + Number(simulationTime,3) + " s / " + (paused ? "PAUSED" : "120 Hz steps"));
    Text(PanelX+22,397,"Grid: 64 x 64     domain length: 4");
    Text(PanelX+22,425,"Viscosity: " + Scientific(viscosity));
    Text(PanelX+22,453,"Dye diffusion: " + Scientific(diffusion));
    Line(PanelX+22,475,PanelX+342,475,{0.16,0.21,0.27});
    const Vec2 normalized = (selected + Vec2{2,2}) / 4.0;
    Text(PanelX+22,504,"PROBE: ("+Number(selected.x,2)+", "+Number(selected.y,2)+")",Muted);
    Text(PanelX+22,533,"Dye: " + Number(fluid.SampleDye(normalized),4),Accent);
    Text(PanelX+22,561,"Speed: " + Number(fluid.SampleVelocity(normalized).Length(),4));
    const double before = stage == 4 ? fluid.DivergenceRms() : fluid.divergenceBeforeProjection;
    const double after = stage == 4 ? fluid.DivergenceRms() : fluid.divergenceAfterProjection;
    Text(PanelX+22,601,stage == 4 ? "Velocity divergence RMS" : "Last projection / divergence RMS",Muted);
    Text(PanelX+22,629,"Before: " + Scientific(before));
    Text(PanelX+22,657,"After:    " + Scientific(after),Accent);
    const double reduction = before > 1e-12 ? 100.0 * (1.0-after/before) : 0.0;
    Text(PanelX+22,685,stage == 4 ? "Pressure solve: not active" : "Reduction: " + Number(reduction,2) + "%",Muted);
    Line(PanelX+22,706,PanelX+342,706,{0.16,0.21,0.27});
    Text(PanelX+22,737,stage == 4 ? "1-4: velocity   [ / ]: diffusion" : "Drag: dye + force   P: projection",Muted);
    Text(PanelX+22,762,"Space: pause   R: reset   F12: capture",Muted);
    Text(PanelX+22,786,stage == 6 ? "B: buoyancy   E: emitters   Tab: next" : "Tab: next stage   +/-: source strength",Muted);
}

void DrawPanel() {
    if (stage >= 4) { DrawFluidPanel(); return; }
    Rect(PanelX, 110, 364, 690, Surface);
    Text(PanelX + 22, 143, stage == 2 ? "GRAVITY SOURCES" : "FIELD PRESETS", Muted);
    const FieldType types[] = {FieldType::Uniform, FieldType::Rotation,
                               FieldType::Source, FieldType::Saddle};
    for (int i = 0; i < 4; ++i) {
        const double y = 158 + i * 47.0;
        if (stage == 2) {
            if (i < static_cast<int>(gravitySources.size())) {
                const auto& source = gravitySources[i];
                Rect(PanelX + 16, y, 332, 38, i == selectedSource ? Color{0.11,0.24,0.25} : Color{0.086,0.12,0.16});
                Text(PanelX + 29, y + 24, std::to_string(i + 1) + "   mass " + Number(source.mass,1)
                     + "   (" + Number(source.position.x,1) + "," + Number(source.position.y,1) + ")");
            }
            continue;
        }
        const bool active = types[i] == fieldType;
        Rect(PanelX + 16, y, 332, 38, active ? Color{0.11, 0.24, 0.25} : Color{0.086, 0.12, 0.16});
        if (active) Rect(PanelX + 16, y, 3, 38, Accent);
        Text(PanelX + 29, y + 24, std::to_string(i + 1) + "   " + FieldName(types[i]),
             active ? Accent : TextColor, GLUT_BITMAP_HELVETICA_18);
    }
    if (stage == 2) {
        Text(PanelX + 22, 361, "a = sum G m r / (r^2 + eps^2)^1.5", Accent);
        Text(PanelX + 22, 384, "G = " + Number(0.4 * strength,2) + "    eps = 0.12    sources: " + std::to_string(gravitySources.size()), Accent);
    } else Text(PanelX + 22, 384, Formula(fieldType), Accent, GLUT_BITMAP_HELVETICA_18);
    Text(PanelX + 22, 414, "Strength k                 " + Number(strength, 2));
    Text(PanelX + 22, 440, stage == 3 ? "Particles                    " + std::to_string(tracers.size()) : "Grid samples              " + std::to_string(gridSize) + " x " + std::to_string(gridSize));
    Text(PanelX + 22, 466, stage == 3 ? "t = " + Number(simulationTime,2) + " s    " + (paused ? "PAUSED" : "RK2 / 120 Hz") : (normalizedArrows ? "Arrow length               normalized" : "Arrow length               compressed magnitude"));
    Line(PanelX + 22, 488, PanelX + 342, 488, {0.16, 0.21, 0.27});
    Text(PanelX + 22, 516, "PROBE / CLICK IN THE FIELD", Muted);
    const Vec2 value = CurrentField(selected);
    const Differential differential = CurrentDerivatives(selected);
    Text(PanelX + 22, 545, "p = (" + Number(selected.x) + ", " + Number(selected.y) + ")");
    Text(PanelX + 22, 571, "F = (" + Number(value.x) + ", " + Number(value.y) + ")", Accent);
    Text(PanelX + 22, 597, "|F| = " + Number(value.Length()));
    Text(PanelX + 22, 633, "Divergence       " + Number(differential.divergence));
    Text(PanelX + 22, 659, "Curl (z)             " + Number(differential.curl));
    Text(PanelX + 22, 687, "Central differences, h = 0.0001", Muted);
    Line(PanelX + 22, 708, PanelX + 342, 708, {0.16, 0.21, 0.27});
    Text(PanelX + 22, 738, stage == 2 ? "Click: add / select   Drag: move   M: mass" : "1-4  preset     +/-  strength     G  grid", Muted);
    Text(PanelX + 22, 762, stage == 2 ? "Right click: remove    R: reset    +/-: G" : "N  normalize     R  reset     F12  capture", Muted);
    Text(PanelX + 22, 786, "F1-F7: stage   Space: pause   F12: capture", Muted);
}

std::wstring CaptureName() {
    if (!captureNameOverride.empty()) return captureNameOverride;
    if (stage == 2) return L"01_\uc911\ub825\ud569\uc131.png";
    if (stage == 3) return L"01_\uc785\uc790\uc774\ub3d9.png";
    if (stage == 4) return L"01_\uc774\ub958.png";
    if (stage == 5) return L"01_\ud22c\uc601\uc801\uc6a9.png";
    if (stage == 6) return L"01_\uc5f0\uae30\ubd84\ucd9c.png";
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
    const char* headings[] = {"", "01    Vector fields", "02    Gravity field", "03    Particle advection", "04    Advection and diffusion", "05    Incompressible fluid", "06    Interactive smoke"};
    Text(54, 77, headings[stage], TextColor, GLUT_BITMAP_HELVETICA_18);
    Text(54, 102, stage >= 4 ? "Semi-Lagrangian transport / implicit diffusion / staggered velocities" : stage == 2 ? "Softened inverse-square acceleration / source superposition" : stage == 3 ? "dx/dt = u(x) / midpoint RK2 / fixed 1/120 second steps" : "Position -> vector -> magnitude -> divergence & curl", Muted);
    Text(PanelX + 165, 45, "STAGE 0" + std::to_string(stage) + " / 07", Muted);
    DrawPlot();
    DrawPanel();
    Text(54, 842, stage == 4 ? "Prescribed velocity transports dye. Pressure projection is introduced in stage 5." : stage >= 5 ? "Dye is a transported tracer.  2D incompressible flow; gas compressibility is not modeled." : stage == 3 ? "Cyan: moving tracers and last 80 substeps.  Arrow field supplies velocity, not force." : "Color: low magnitude (cyan) -> high (amber).  Display length does not change F.", Muted);
    Text(PanelX, 829, captureStatus, Muted);
    if (captureRequested) {
        captureRequested = false;
        const bool saved = SaveFrame();
        captureStatus = saved ? "PNG saved. F12 replaces this preset capture." : "Capture failed. Check the output folder.";
        std::printf("Capture %s (stage %d, t=%.3f, k=%.2f)\n", saved ? "saved" : "FAILED", stage, simulationTime, strength);
        std::fflush(stdout);
        if (quitAfterCapture) { glutSwapBuffers(); std::exit(saved ? 0 : 2); }
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
    case 'm': case 'M':
        if (stage == 2 && selectedSource >= 0 && selectedSource < static_cast<int>(gravitySources.size())) {
            double& mass = gravitySources[selectedSource].mass;
            mass = mass >= 4.0 ? 0.5 : mass + 0.5;
        }
        break;
    case 'p': case 'P': if (stage >= 5) projectionEnabled = !projectionEnabled; break;
    case 'b': case 'B': if (stage == 6) buoyancyEnabled = !buoyancyEnabled; break;
    case 'e': case 'E': if (stage == 6) emitterEnabled = !emitterEnabled; break;
    case '[': diffusion = std::max(0.0,diffusion * 0.5); break;
    case ']': diffusion = std::min(0.1,std::max(0.00001,diffusion * 2.0)); break;
    case ' ': paused = !paused; break;
    case '\t': stage = stage == AvailableStages ? 1 : stage + 1; ResetStage(); break;
    case 'r': case 'R': ResetStage(); break;
    }
    if (stage >= 4 && key >= '1' && key <= '4') fluid.SetPresetVelocity(fieldType,strength);
    if (stage == 4 && (key == '+' || key == '=' || key == '-' || key == '_')) fluid.SetPresetVelocity(fieldType,strength);
    glutPostRedisplay();
}

void Special(int key, int, int) {
    if (key >= GLUT_KEY_F1 && key < GLUT_KEY_F1 + AvailableStages) { stage = key - GLUT_KEY_F1 + 1; ResetStage(); glutPostRedisplay(); }
    if (key == GLUT_KEY_F12) { captureRequested = true; glutPostRedisplay(); }
}

void Mouse(int button, int state, int x, int y) {
    if (state == GLUT_UP) { draggingSource = false; mouseHeld = false; return; }
    const double px = (x - viewportX) / viewportScale;
    const double py = (y - viewportYTop) / viewportScale;
    if (px >= PlotX && px <= PlotX + PlotSize && py >= PlotY && py <= PlotY + PlotSize) {
        selected = {(px - PlotX) / PlotSize * 4.0 - 2.0,
                     2.0 - (py - PlotY) / PlotSize * 4.0};
        if (stage >= 4) {
            mouseHeld = button == GLUT_LEFT_BUTTON;
            lastMousePosition = selected;
            fluid.AddDye((selected.x+2)/4.0,(selected.y+2)/4.0,0.8,0.05);
        }
        if (stage == 2) {
            int nearest = -1;
            double distance = 0.24;
            for (int i = 0; i < static_cast<int>(gravitySources.size()); ++i) {
                const double d = (selected - gravitySources[i].position).Length();
                if (d < distance) { distance = d; nearest = i; }
            }
            if (button == GLUT_RIGHT_BUTTON && nearest >= 0) {
                gravitySources.erase(gravitySources.begin() + nearest);
                selectedSource = gravitySources.empty() ? -1 : 0;
            } else if (button == GLUT_LEFT_BUTTON) {
                if (nearest < 0 && gravitySources.size() < 16) {
                    gravitySources.push_back({selected, 1.0}); nearest = static_cast<int>(gravitySources.size()) - 1;
                }
                selectedSource = nearest; draggingSource = nearest >= 0;
            }
        }
    } else if (stage < 4 && px >= PanelX + 16 && px <= PanelX + 348) {
        for (int i = 0; i < 4; ++i) {
            if (py >= 158 + i * 47.0 && py <= 196 + i * 47.0) {
                if (stage == 2 && i < static_cast<int>(gravitySources.size())) selectedSource = i;
                else if (stage != 2) fieldType = static_cast<FieldType>(i);
            }
        }
    }
    glutPostRedisplay();
}

void Motion(int x, int y) {
    const double px = (x - viewportX) / viewportScale;
    const double py = (y - viewportYTop) / viewportScale;
    Vec2 p{std::max(-2.0,std::min(2.0,(px-PlotX)/PlotSize*4.0-2.0)),
           std::max(-2.0,std::min(2.0,2.0-(py-PlotY)/PlotSize*4.0))};
    if (stage >= 4 && mouseHeld) {
        const Vec2 delta = p - lastMousePosition;
        if (stage >= 5) fluid.AddVelocity((p.x+2)/4.0,(p.y+2)/4.0,delta.x*5.0*strength,delta.y*5.0*strength,0.06);
        lastMousePosition = p; selected = p; glutPostRedisplay(); return;
    }
    if (stage != 2 || !draggingSource || selectedSource < 0) return;
    gravitySources[selectedSource].position = p;
    selected = p;
    glutPostRedisplay();
}

void Timer(int) {
    const int now = glutGet(GLUT_ELAPSED_TIME);
    const double elapsed = std::max(0, now - lastTimerTime) / 1000.0;
    recordElapsed += elapsed;
    lastTimerTime = now;
    if (stage >= 3 && !paused) {
        accumulator += std::min(0.1,elapsed);
        const double dt = 1.0/120.0;
        while (accumulator >= dt) {
            if (stage == 3) for (auto& tracer : tracers) {
                tracer.position = AdvectRK2(tracer.position, dt, fieldType, strength);
                bool wrapped = false;
                if (tracer.position.x > 2) { tracer.position.x -= 4; wrapped = true; }
                if (tracer.position.x < -2) { tracer.position.x += 4; wrapped = true; }
                if (tracer.position.y > 2) { tracer.position.y -= 4; wrapped = true; }
                if (tracer.position.y < -2) { tracer.position.y += 4; wrapped = true; }
                if (wrapped) tracer.trail.clear();
                tracer.trail.push_back(tracer.position);
                if (tracer.trail.size() > 80) tracer.trail.pop_front();
            }
            if (stage >= 4) {
                if (mouseHeld) fluid.AddDye((selected.x+2)/4.0,(selected.y+2)/4.0,24.0*dt,0.05);
                if (stage == 4) fluid.StepDye(dt,diffusion);
                else {
                    if (stage == 6 && emitterEnabled) {
                        fluid.AddDye(0.47+0.05*std::sin(0.7*simulationTime),0.12,8.0*dt,0.055);
                        fluid.AddVelocity(0.47+0.05*std::sin(0.7*simulationTime),0.12,0.6*std::sin(1.3*simulationTime)*dt,3.0*dt,0.06);
                        fluid.AddDye(0.68,0.14,5.0*dt,0.055);
                        fluid.AddVelocity(0.68,0.14,-0.7*std::sin(0.9*simulationTime)*dt,2.4*dt,0.06);
                    }
                    fluid.Step(dt,viscosity,diffusion,projectionEnabled,true,stage == 6 && buoyancyEnabled);
                }
            }
            simulationTime += dt; accumulator -= dt;
        }
        glutPostRedisplay();
    }
    if (captureAfter >= 0.0 && (stage >= 3 ? simulationTime : recordElapsed) >= captureAfter) {
        captureAfter = -1.0; captureRequested = true; glutPostRedisplay();
    }
    glutTimerFunc(16, Timer, 0);
}
} // namespace

int Run(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--self-test") == 0) return RunFieldMathTests() | RunDynamicsTests() | RunFluidGridTests();
        if (std::strcmp(argv[i], "--legacy") == 0) { LegacyMain(argc, argv); return 0; }
    }
    int wideCount = 0;
    LPWSTR* wideArgs = CommandLineToArgvW(GetCommandLineW(), &wideCount);
    if (wideArgs) {
        for (int i = 1; i + 1 < wideCount; ++i) {
            if (std::wstring(wideArgs[i]) == L"--capture-dir") captureDirectory = wideArgs[++i];
            else if (std::wstring(wideArgs[i]) == L"--stage") stage = std::max(1,std::min(AvailableStages,_wtoi(wideArgs[++i])));
            else if (std::wstring(wideArgs[i]) == L"--capture-after") captureAfter = _wtof(wideArgs[++i]);
            else if (std::wstring(wideArgs[i]) == L"--capture-name") captureNameOverride = wideArgs[++i];
            else if (std::wstring(wideArgs[i]) == L"--strength") requestedStrength = _wtof(wideArgs[++i]);
        }
        for (int i = 1; i < wideCount; ++i) {
            if (std::wstring(wideArgs[i]) == L"--quit-after-capture") quitAfterCapture = true;
            if (std::wstring(wideArgs[i]) == L"--no-projection") requestedNoProjection = true;
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
    ResetStage();
    if (requestedStrength > 0.0) strength = std::max(0.25,std::min(3.0,requestedStrength));
    if (stage == 4 && requestedStrength > 0.0) fluid.SetPresetVelocity(fieldType,strength);
    if (requestedNoProjection) projectionEnabled = false;
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    glutKeyboardFunc(Keyboard);
    glutSpecialFunc(Special);
    glutMouseFunc(Mouse);
    glutMotionFunc(Motion);
    lastTimerTime = glutGet(GLUT_ELAPSED_TIME);
    glutTimerFunc(16, Timer, 0);
    // Event-driven drawing avoids an idle loop consuming a CPU core.
    glutMainLoop();
    return 0;
}
} // namespace fieldlab

int main(int argc, char** argv) { return fieldlab::Run(argc, argv); }
