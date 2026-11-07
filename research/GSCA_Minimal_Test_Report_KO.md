# GSCA: auxiliary replacement와 outcome indicators만 변경한 테스트

요청한 두 변경만 적용하고 원래 GSCA 코드를 다시 실행했다. 기존 ALS 추정, lambda 선택, single bootstrap, parameter 신뢰구간 및 percentile 계산은 유지했다. **Nested bootstrap과 variance correction은 포함하지 않았다.**

## 변경한 계산

1. 각 bootstrap의 in-bag component scores에서 focal predictor를 같은 outcome 방정식의 나머지 component predictors에 회귀시킨다. 절편을 포함한 ordinary least squares를 사용한다. 이 식을 OOB scores에 적용한 예측값으로 focal score를 대체한다.
2. B 경로의 예측 outcome은 target component에 연결된 관측 indicators로 바꾼다. 부모 component scores → B → C로 indicators를 예측하고, in-bag 표준화를 역변환하여 원래 관측 단위에서 오차를 계산한다. Indicator별 제곱오차를 먼저 구한 뒤 indicator와 OOB 관측치에 걸쳐 평균한다.

각 bootstrap에서 delta OPE는 replacement prediction의 평균제곱오차에서 full prediction의 평균제곱오차를 뺀 값이다. 기존과 같이 bootstrap별 결과를 동일하게 가중하여 평균한다. W/C/B를 다시 적합하지 않고 작은 auxiliary regression만 추가하므로 GSCA 적합 횟수는 bootstrap 수에 비례한다.

관측 outcome 자체가 추정된 target weights에 따라 달라지는 문제를 없앤 것이다. 예측치에는 학습된 W/C/B의 영향이 계속 반영된다. 현재 구현 범위는 서로 겹치지 않는 indicator blocks이며, 각 방정식의 parent blocks를 예측 시점에 관측할 수 있다고 가정한다. C 경로는 기존 standardized-indicator reconstruction loss를 유지한다.

## 원본 보존 검증: 통과

Native MATLAB R2026a에서 원본 main commit `9d586446970a47ebd420516a0451143ba608fff9`와 수정본을 같은 seed로 비교했다.

- 원본 INI의 모든 기존 필드, W/C/B, CV 결과와 선택된 lambda, 기존 GoF가 정확히 일치했다.
- 모든 bootstrap parameter draws, parameter estimate/SE/CI, 최종 난수 상태가 정확히 일치했다.
- 원래 percentile 순위와 nonpositive-draw 비율 계산을 유지했다.
- Auxiliary OLS 산술, raw-unit indicator 예측, OOB target 입력 제외, sign invariance, 단일 indicator OLS에서의 FWL 관계를 검증했다.
- 기존 ACSI replica도 같은 원본 solver로 정상 실행했다.

원본과의 동일성 비교는 serial 실행에서 검증했다. 코드의 serial/parfor 두 경로에 같은 evaluation helper를 연결했으며, 이번 simulation은 serial 방식으로 실행했다.

## 작은 simulation

P=3, 각 block의 J_p=3. 두 predictor components의 상관은 .30이고, outcome의 다른 predictor 계수는 .35로 고정했다. Within-block indicator 상관은 .50이며 indicator의 모집단 분산은 1이다. 정규분포 exact-composite DGP를 사용했다.

N=80/200 × lambda=0/.1 × focal beta=0/.1/.3의 12개 조건에서 각 100개 데이터셋, 각 데이터셋당 199회 bootstrap을 실행했다. **총 1,200개 데이터셋 모두 완료했고 실패는 0개였다.** Bootstrap GSCA 적합은 238,800회이며 initial fit은 별도이다. 조건별 데이터 생성 seed와 bootstrap seed를 모두 저장했다.

기각은 원본 TABLE.B의 column 5, 즉 one-sided percentile lower bound가 0보다 큰 경우로 계산했다. B=199에서 원본의 rounded order-statistic 규칙을 그대로 사용했다.

| N | lambda | beta | 평균 delta OPE | 평균의 MCSE | 평균 indicator OPE | 기각률 | 기각률의 95% Wilson 구간 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 80 | 0 | 0 | -0.02005 | 0.00100 | 0.98228 | 0% | 0.0–3.7% |
| 80 | 0 | 0.1 | -0.01122 | 0.00179 | 0.95442 | 0% | 0.0–3.7% |
| 80 | 0 | 0.3 | 0.04217 | 0.00365 | 0.86327 | 1% | 0.2–5.4% |
| 80 | 0.1 | 0 | -0.01313 | 0.00156 | 0.95070 | 0% | 0.0–3.7% |
| 80 | 0.1 | 0.1 | -0.00235 | 0.00245 | 0.93985 | 0% | 0.0–3.7% |
| 80 | 0.1 | 0.3 | 0.04187 | 0.00439 | 0.85669 | 4% | 1.6–9.8% |
| 200 | 0 | 0 | -0.00681 | 0.00033 | 0.93656 | 0% | 0.0–3.7% |
| 200 | 0 | 0.1 | 0.00021 | 0.00093 | 0.91545 | 0% | 0.0–3.7% |
| 200 | 0 | 0.3 | 0.05073 | 0.00271 | 0.83032 | 15% | 9.3–23.3% |
| 200 | 0.1 | 0 | -0.00488 | 0.00039 | 0.94350 | 0% | 0.0–3.7% |
| 200 | 0.1 | 0.1 | 0.00298 | 0.00110 | 0.90592 | 0% | 0.0–3.7% |
| 200 | 0.1 | 0.3 | 0.04864 | 0.00239 | 0.82427 | 36% | 27.3–45.8% |

### 이번 결과에서 확인한 점

- beta=0에서는 네 조건 모두 평균 delta OPE가 음수였다.
- beta=.1에서는 N=80의 평균 delta OPE가 음수였고, N=200에서는 약 0.00021과 0.00298이었다.
- beta=.3에서는 평균 delta OPE가 네 조건 모두 양수였다. 하지만 기존 percentile 방식의 기각률은 N=80에서 1%/4%, N=200에서 15%/36%였다(lambda=0/.1 순서).
- 따라서 두 변경 후의 계산은 정상 작동했지만, 평균 delta OPE가 양수인 조건에서도 기존 검정의 기각률이 높지는 않았다. 이 결과만으로 검정의 calibration이나 estimator bias가 개선되었다고 결론 내리지는 않는다.

여기서 beta는 데이터 생성 계수이다. 이번 테스트는 finite-N predictive-null 경계를 별도로 구하지 않았으므로 위 비율은 조건별 **기각률**로 제시했다. 100회 중 기각 0회는 모집단 기각확률이 정확히 0이라는 뜻이 아니다. MCSE와 Wilson 구간은 simulation 불확실성을 나타낸다.

## ACSI replica

원본 예제의 W0/C0/B0를 유지하고 lambda=0, bootstrap=199로 실행했다. 아래 delta와 lower bound는 raw indicator units의 평균제곱오차 차이이다. Outcome block 간 단위가 다르면 효과 크기의 수치도 직접 비교하기 어렵다.

| 경로 | 경로계수 | 평균 delta OPE | one-sided lower bound | nonpositive-draw 비율 |
|---|---:|---:|---:|---:|
| 1 → 2 | 0.9177 | 416.847 | 149.130 | 0.0000 |
| 1 → 3 | -0.0874 | -3.401 | -14.687 | 0.6231 |
| 2 → 3 | 0.9570 | 55.798 | 9.189 | 0.0201 |
| 1 → 4 | -0.0399 | -0.490 | -2.209 | 0.6935 |
| 2 → 4 | 0.7321 | 22.380 | 7.773 | 0.0000 |
| 3 → 4 | 0.3145 | 9.125 | 3.679 | 0.0000 |
| 4 → 5 | -0.4697 | 89.764 | -61.437 | 0.1206 |
| 4 → 6 | 0.4943 | 64.191 | 22.893 | 0.0151 |
| 5 → 6 | -0.4410 | 55.909 | 19.715 | 0.0050 |

기존 component-loss 평균도 별도로 저장했다. Outcome과 단위가 달라졌으므로 기존 평균과 새로운 indicator 평균의 단순 차이를 bias correction 크기로 해석하지 않는다.

## 코드와 재현

`GSCA_Reg.m`은 기존 인수를 그대로 받으며, 같은 경로에 `GSCA_Indicator_OOB.m`을 두면 된다. 새 통계량은 `ETC.delOPE_b_Boot`, `ETC.indicatorOPE_Boot`, `INI.indicatorOPE_byComponent`에서 확인할 수 있다. 기존 component 기준 draws는 `ETC.delOPE_b_legacy_Boot`에 보관했다. 기존 `INI.GoF`는 이전 의미를 유지한다.

```matlab
addpath('research');
test_gsca_minimal_update;
test_gsca_original_flow;
Run_Minimal_Indicator_Test(80);
Run_Minimal_Indicator_Test(200);
run_minimal_indicator_example;
```

MATLAB과 원본 solver의 zscore용 Statistics and Machine Learning Toolbox가 필요하다.

- Numerical commit: `0ef5a104e62bd7ad52164c0c7580e5117216301d`
- [완료된 실행과 검증 로그](https://github.com/PsycheMatrica/GSCA.Reg_Prime/actions/runs/42339690205)
- [최소 변경 PR #2](https://github.com/PsycheMatrica/GSCA.Reg_Prime/pull/2)
- 연구 branch: `research/indicator-auxiliary-minimal-20261106`

결과 패키지는 수정 코드, 원본 solver, 모든 replicate 결과와 seed, MAT/CSV/JSON, ACSI 결과 및 파일 hashes를 포함한다.
