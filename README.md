# VectorField

C++ / OpenGL로 벡터장에서 시작해 2차원 유체 시뮬레이션까지 단계별로 개발하는 프로젝트입니다.

## 단계별 결과

| 단계 | 내용 | 상태 |
| --- | --- | --- |
| [1.벡터장](1.벡터장/README.md) | 4가지 장, 선택 위치의 벡터·발산·회전, 수학 검증, 실제 화면 캡처 | 완료 |
| 2.중력장 | 여러 중력원의 역제곱 중력과 중첩 | 다음 단계 |
| 3.입자이동 | RK2로 속도장을 따라 움직이는 추적 입자 | 예정 |
| 4.이류와확산 | 격자의 염료 이류·확산 | 예정 |
| 5.비압축유체 | 압력 투영을 포함한 2D 유체 솔버 | 예정 |
| 6.상호작용유체 | 힘과 연기 주입, 부력, 실시간 상호작용 | 예정 |

각 번호 폴더에는 개발 설명, 계산과정, 검증결과, 실제 캡처와 해당 시점의 실행본·소스를 보관합니다. 전체 기준은 [단계별 개발계획](단계별_개발계획.md)에 정리했습니다. 현재 유체 시뮬레이션은 구현 전이며 1단계부터 순서대로 진행합니다.

## 실행

`run.cmd`를 더블 클릭하세요. 기본 실행 파일은 `x64\Release\opgnel_project_smoke.exe`입니다. 실행 파일이 없으면 먼저 빌드합니다.

기본 실행은 1단계 벡터장입니다. `1`~`4`로 장을 선택하고, 마우스로 위치를 선택해 값을 확인합니다. `+`/`-`는 강도, `G`는 격자 수, `N`은 길이 정규화, `R`은 초기화, `F12`는 PNG 저장입니다. 미분은 격자와 별도로 중앙차분으로 계산합니다.

단계별 보관본은 `1.벡터장/실행.cmd`로 실행합니다. 원래 저장소의 `Smoke` 중력장 데모는 `run.cmd --legacy`로 열 수 있습니다. 기존 데모의 수학적 수정은 2단계에서 진행합니다.

## 빌드

다음 도구가 필요합니다.

- Git for Windows
- Visual Studio 2022의 **C++를 사용한 데스크톱 개발** 워크로드
- MSVC v143, Windows 10 또는 11 SDK, **Windows용 C++ CMake 도구** 구성 요소

프로젝트 폴더에서 실행하세요.

```powershell
.\build.cmd
.\run.cmd
```

`build.cmd`는 Visual Studio 설치 위치를 자동으로 찾고, 공식 freeglut 3.6.0을 `.deps` 폴더에 내려받아 빌드한 뒤 앱을 빌드합니다. 최초 빌드에는 인터넷 연결이 필요합니다. freeglut 소스는 커밋 `566cd459201ae29003ef5686c06506d8646e55c6`으로 고정되어 있습니다.

기본 구성은 `Release / x64`입니다. 다른 구성은 다음처럼 지정합니다.

```powershell
.\build.cmd -Configuration Debug
.\build.cmd -Platform Win32
```

x64 실행 파일은 `x64\Release` 또는 `x64\Debug`에, Win32 실행 파일은 `Release` 또는 `Debug`에 생성됩니다. 필요한 `freeglut.dll`은 실행 파일 옆에 자동 복사됩니다. `run.cmd`는 기본 `Release / x64` 실행 파일을 실행합니다.

## Visual Studio 2022에서 열기

먼저 사용할 플랫폼으로 `build.cmd`를 한 번 실행해 freeglut를 준비한 뒤 `opgnel_project_smoke.sln`을 여세요. 구성은 `Release` 또는 `Debug`, 솔루션 플랫폼은 `x64` 또는 `x86`을 선택하고 빌드·실행하면 됩니다. 솔루션의 `x86`은 프로젝트의 `Win32`에 대응합니다.

현재 진입점과 화면은 `opgnel_project_smoke/FieldLab.cpp`, 벡터장 계산은 `FieldMath.h`에 있습니다. 원래 `source.cpp`는 `LegacyMain()`으로 보존했습니다. `x64/Release/opgnel_project_smoke.exe --self-test`로 수학 검증을 실행합니다.
