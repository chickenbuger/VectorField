# 결합 연산 독립 검증

이 검증기는 `TwoMaterialTests.h`와 별도로 작성했습니다. `WaterMomentum()`, `GasMomentum()`, `TotalKineticEnergy()` 같은 진단 합계를 사용하지 않고, 모든 내부 MAC 면의 실제 속도와 모든 물 입자의 속도를 읽어 운동량과 운동에너지를 다시 합산합니다.

## 검사한 범위

다음 조합의 64개 상태를 고정 난수 시드로 만들었습니다. 물 입자는 내부뿐 아니라 네 벽과 네 모서리 가까이에도 배치합니다.

| 변수 | 값 |
| --- | --- |
| 격자 한 축 해상도 | 8, 12, 24, 48 |
| 기체 기준 밀도 | 0.1, 1.2, 80, 1000 |
| 결합 연산 시간 간격 | 0.00001, 1/120, 0.2, 10초 |

각 상태에서 `ExchangeOnly(dt)` 전후의 두 매체 운동량 합과 운동에너지를 비교합니다. 0.2초와 10초는 이 결합 연산의 강한 감쇠에서도 오버슈트가 없는지 확인하는 입력이며, 전체 유체 시간 적분에 사용할 권장 간격이 아닙니다.

추가로 공동 속도 평형과 Off·OneWay·TwoWay의 상태 변화도 확인합니다. Off에서는 두 속도가 유지되고, OneWay에서는 기체 속도가 유지되며 물만 바뀌고, TwoWay에서는 반작용까지 적용되어야 합니다.

실제 실행 결과:

```text
Independent review scenarios: 64 randomized near-wall cases + equilibrium + 3 modes
Maximum relative actual-state momentum residual: 1.4551874580857375e-15
Maximum relative kinetic energy increase: 0
Independent exchange review failures: 0
```

상대 운동량 잔차는 `|ΔP물+ΔP기체| / (1+|P물,전|+|P기체,전|)`입니다. 에너지 증가는 `(E후−E전)/(1+E전)`의 양수 부분 중 최대값입니다. 증가 0은 에너지가 일정하다는 뜻이 아니라, 검사한 모든 상태에서 에너지가 증가하지 않았다는 뜻입니다.

이 검사는 **결합 교환 연산만** 검증합니다. 전체 `Step()`에는 중력, 분출원, 압력 투영, 물의 위치·벽 보정이 있어 전체 장면의 운동량은 별도로 바뀔 수 있습니다. 연기 농도 제외는 시각화용 추적량 처리이며, 불침투 경계나 기체 질량 보존의 검증은 아닙니다.

## 기본 해상도 전체 진행 확인

UI와 같은 기본 48×48 격자, 400개 물 입자, 기체 기준 밀도 80, 분출원 켬, TwoWay 설정도 별도로 실행했습니다. 600회 고정 시간 단계로 5초를 계산했으며 모든 단계에서 유한값, 물 입자의 벽 안 유지, 기체 벽의 법선 속도 0, 물 입자 수·질량 유지를 확인했습니다.

| 측정 | 실제 결과 |
| --- | --- |
| 600스텝 실행 시간 | 약7.664초 — 실행한 PC의 Release 빌드 실측 |
| 각 교환 연산의 최대 절대 운동량 잔차 | 5.2618524848329024×10⁻¹² |
| 각 교환 연산의 최대 상대 잔차 | 3.2528772826823303×10⁻¹³ |
| 1초의 Off/TwoWay 물 입자 위치 RMS 차이 | 0.34329252499077401 |
| 같은 조건의 내부 MAC 면 기체 속도 RMS 차이 | 0.67849435436651406 |

전체 진행의 상대 잔차는 `|ΔP물+ΔP기체| / (1+|ΔP물|+|ΔP기체|)`로 정규화한 교환 진단값입니다. 앞의 64개 검사에서 사용하는 초기 운동량 기준 정규화와는 다릅니다. Off/TwoWay의 위치와 기체 속도 차이는 결합을 켰을 때 두 매체의 계산 상태가 모두 달라지는지 확인하는 비교입니다. 원인이 분출원만인지, 움직이는 물과 정지 기체 사이 항력까지 포함하는지는 별도의 분출원 켬/끔 대조군으로 구분해야 합니다.

전체 실제 출력은 상위 폴더의 [독립검증결과.txt](../독립검증결과.txt)에 보관했습니다. 이 추가 검사도 종료 코드 0이며, 앞서 말한 보존 범위는 여전히 교환 연산에 한정합니다.

## 다시 실행하기

저장소 루트의 PowerShell에서 Visual Studio 2022에 포함된 CMake를 사용합니다.

```powershell
$vsRoot = & "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products '*' -version '[17.0,18.0)' -property installationPath
$cmake = Join-Path $vsRoot 'Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe'
& $cmake -S '8.다중재료상호작용/검증' -B build/stage8-independent-review -G 'Visual Studio 17 2022' -A x64
& $cmake --build build/stage8-independent-review --config Release
& './build/stage8-independent-review/Release/independent_coupling_review.exe'
```

사용한 표준은 C++14이며 OpenGL이나 freeglut는 필요하지 않습니다. 종료 코드 0은 위 검증 범위를 통과했다는 뜻입니다.
