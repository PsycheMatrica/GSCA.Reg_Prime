# Indicator 기반 GSCA 보정: 점추정 결과와 남은 문제

## 현재 결론

GSCA 구조 경로의 outcome을 추정된 성분에서 **고정된 관측지표 블록**으로 바꾸고, 지표별 제곱오차의 평균을 평가하도록 구현했다. MATLAB R2026a에서 독립 검증을 통과한 뒤 기본 시뮬레이션 96개 조건/reference 행, 원자료 수준 분석 9,600회를 완료했다. 동일 원자료의 oracle/sample-reference 비교가 포함되므로 9,600개의 독립 실험이라는 뜻은 아니다.

보정 후 관측된 절대 편향은 96개 중 89개 행에서, RMSE는 90개 행에서 줄었다. 다중 지표($J_p=3$) 조건만 보면 각각 48개 중 45개 행에서 줄었다. 이는 관측된 simulation 성능이며, 불편성이나 일양적인 우월성의 증명이 아니다.

**Nested bootstrap은 현재의 해결책으로 채택하지 않는다.** 바깥 재표집은 최종 추정량의 변동을 추정하려는 절차지만, bootstrap/OOB가 평가하는 학습 위험과 원래 finite-$N$ target의 차이를 자동으로 보정하지 않는다. 진행 중이던 다섯 nested 실험은 중단했으며, 그 부분 결과를 검증 근거로 사용하지 않았다. 새 indicator 추정량의 기본 설정도 Outer=0이다.

## 공통 target과 예측 규칙

연구의 연결은 MLR under MVN → bootstrap MLR 및 single-outcome ERA → predictor/outcome 블록을 갖는 GSCA이다. 모든 단계에서 평가 outcome은 고정된 관측변수이며, components는 예측 알고리즘 안에서 학습되는 표현이다.

목표 성분 $q$에 속한 지표가 $J_q$개라면

$$
\mathrm{EPE}^{(q)}_m(N)=\frac1{J_q}
E_{\mathcal D_N,Z}\sum_{j=1}^{J_q}
\{Y_{qj}-\widehat f_{m,qj}(X;\mathcal D_N)\}^2,
\qquad
\Delta\mathrm{EPE}^{(q)}_p(N)=
\mathrm{EPE}^{(q)}_{-p}(N)-\mathrm{EPE}^{(q)}_1(N).
$$

부모 블록의 score에 구조계수 $B$를 적용해 목표 성분을 **예측**하고, loading $C$와 훈련표본의 척도 변환을 통해 원래 단위의 지표 예측값을 만든다. 실제 목표 지표로 계산한 성분 score를 자기 지표 예측에 사용하지 않는다. 각 목표 블록의 부모 지표가 예측 시점에 관찰된다는 local conditional prediction이다.

Predictor replacement는 해당 source score를 같은 구조식의 다른 부모 score들에 대한 affine linear projection으로 대체한다. Outcome 식은 다시 적합하지 않는다. 참값은 fitted weights가 주어졌을 때의 **모집단 projection**을 사용한다. Sample reference는 그 projection을 전체 원자료의 predictor moments로 추정하는 nuisance approximation이며, 다른 target으로 재정의하지 않는다.

현재 구현은 simple disjoint reflective blocks와 고정된 $\lambda$를 지원한다. Crossloading, 자기 경로, 목표 지표의 predictor 입력 사용은 명시적으로 거부한다. 원척도의 손실을 사용하므로 서로 다른 측정 단위를 섞은 평균이나 서로 다른 outcome 블록 간 효과 크기 비교는 별도의 고정된 평가 척도가 필요하다.

## 추정량과 correction의 정확한 의미

각 관측값에 동일한 총 가중치를 주고, 그 관측값이 OOB인 반복들을 먼저 평균한다. $K_i$는 관측값 $i$의 OOB 평가 횟수다.

$$
\widehat{\Delta\mathrm{EPE}}^{\rm raw}_p
=\frac1{NJ_q}\sum_{i,j}
\overline{e_{-p,bij}^2-e_{1,bij}^2},
\qquad
\widehat C_p=\frac1{NJ_q}\sum_{i,j}
\{s_b^2(e_{1,bij})-s_b^2(e_{-p,bij})\}.
$$

보정 후보는 $\widehat{\Delta\mathrm{EPE}}^{\rm VC}_p
=\widehat{\Delta\mathrm{EPE}}^{\rm raw}_p+\widehat C_p$이다. 분산은 $K_i-1$을 분모로 계산한다. 고정된 outcome에서는 잔차 분산과 예측값 분산이 같으므로 moving-target covariance 항이 사라진다. **지표별 제곱오차와 분산을 계산한 후 평균하며, 평균 잔차를 제곱하지 않는다.**

다음 관계는 정확한 유한-bootstrap 항등식이다.

$$
\widehat{\Delta\mathrm{EPE}}^{\rm VC}_p
=\frac1{NJ_q}\sum_{i,j}
(\bar e_{-p,ij}^{\,2}-\bar e_{1,ij}^{\,2})
+\frac1{NJ_q}\sum_{i,j}
\frac{s_b^2(e_{1,bij})-s_b^2(e_{-p,bij})}{K_i}.
$$

이 항등식은 보정이 bootstrap 예측을 평균한 bagged predictor의 손실 차이와 연결됨을 보여준다. **독립적인 새 훈련표본 $N$개로 적합한 원래 GSCA 위험을 불편 추정한다는 결론은 나오지 않는다.** Bootstrap 훈련표본의 중복 및 OOB 조건화, 비선형 weight 추정과 벌점의 영향을 target과 연결하는 검증이 남아 있다. 바깥 bootstrap을 추가해도 이 연결이 자동으로 성립하지 않는다.

현재 raw/shifted percentile은 진단용 benchmark다. 모두 단측 95% 하한이며, 경험적 $p<.05$와 하한 $>0$을 일치시켰다. 점추정 개선을 근거로 이를 보정된 predictive significance test라고 제시하지 않는다.

## 완료된 기본 실험

$N=80,200$, $J_p=1,3$, $\lambda=0,.1$, 정규/이분산, 계수 0/경계 후보/.30, oracle/sample reference를 교차했다. 각 행은 원자료 100개와 OOB 199회다. OOB 적합 1,910,400회가 완료되었고, 별도 원자료 적합·경계 calibration·참값 계산은 이 수에 포함하지 않았다.

참값은 독립 훈련표본 1,500개의 affine prediction risk를 알려진 모집단 평균·공분산으로 정확히 적분해 평균했다. Test-sample MC 오차는 없지만 training-risk 평균의 MCSE는 남는다. 일반 GSCA 경계는 500개 훈련표본의 common-seed calibration으로 찾은 후보이며, 독립적인 참값 추정으로 검토했다.

정규·무벌점·단일 지표에서는 기존의 정확한 MLR 공식

$$
\Delta\mathrm{EPE}_p(N)=.91b^2-\frac{.8775-b^2-.21b}{N-4}
$$

을 유지한다. 여섯 N/effect 조합의 독립 MC 확인값은 모두 analytic target에서 1.80 MCSE 이내였다.

### 다중 지표·sample reference의 경계 후보 결과

괄호는 target MCSE이다. 양의 참값이 확인된 행의 기각률은 Type I error로 해석하지 않는다.

| N | 자료 | λ | 참값 추정 (MCSE) | Raw bias | VC bias | Raw RMSE | VC RMSE |
|---:|---|---:|---:|---:|---:|---:|---:|
| 80 | 정규 | 0 | 0.000093 (0.000282) | -0.00710 | 0.00154 | 0.02702 | 0.02549 |
| 80 | 이분산 | 0 | -0.001206 (0.000537) | -0.01514 | -0.00137 | 0.03759 | 0.03198 |
| 80 | 정규 | 0.1 | -0.000174 (0.000242) | -0.00558 | 0.00098 | 0.02168 | 0.02101 |
| 80 | 이분산 | 0.1 | 0.000636 (0.000422) | -0.01377 | -0.00271 | 0.03396 | 0.02920 |
| 200 | 정규 | 0 | 0.000035 (0.000109) | -0.00319 | 0.00006 | 0.00817 | 0.00755 |
| 200 | 이분산 | 0 | 0.000942 (0.000190) | -0.00583 | -0.00072 | 0.01410 | 0.01233 |
| 200 | 정규 | 0.1 | 0.000251 (0.000090) | -0.00359 | -0.00112 | 0.00669 | 0.00557 |
| 200 | 이분산 | 0.1 | -0.000108 (0.000167) | -0.00526 | -0.00068 | 0.01443 | 0.01300 |

다중 지표의 경계 후보 16개 reference 행에서 raw와 shifted percentile의 관측 기각률은 모두 0이었다. 특히 일부 후보는 참값이 작은 양수로 확인되었는데도 그러했다. **점추정 보정만으로 개별 OOB percentile 검정의 보수성이 해결된 것은 아니다.**

다중 지표의 성능 비개선 세 행은 모두 $N=200,\lambda=0,b=.30$이었다: 정규/oracle, 이분산/oracle, 이분산/sample. 원자료 100개에서 관찰한 bias/RMSE 차이이며, 소폭 악화를 통계적으로 확정한 결과는 아니다. 전체 행과 MCSE는 요약 파일에 보존했다.

## 구현 검증과 ACSI 적용

직접 계산한 raw-data OLS bootstrap과 GSCA 단일 지표 결과를 비교했으며, affine prediction·projection·지표 평균·부호/단위 변환·OOB 관측값 가중·고정 outcome 분산 항등식 검사를 통과했다. 새 map은 목표 지표가 predictor 입력에 들어가는 구조를 거부한다. 기본 시뮬레이션 9,600회 분석에서 원자료 분석 실패 및 OOB 적합 실패는 모두 0회였다. 관찰된 0회 실패가 모든 모집단에서의 수렴이나 위험의 적분 가능성을 증명하지는 않는다.

기존 ACSI replica($N=150$, 성분 6개, 지표 14개)의 9개 구조 경로에도 $\lambda=0$, OOB 199회로 적용했고 실패가 없었다. 이는 기존 replica에 대한 구현 적용 확인이며 새로운 실증 자료가 아니다. Outer=0이고 percentile은 진단용이므로 검증된 predictive p-value로 해석하지 않는다.

## 다음 검증의 우선순위

1. 단일 bootstrap/OOB 추정량이 원래 finite-$N$ target과 어떻게 다른지 MLR의 analytic truth로 분리해 확인한다.
2. 현재 variance correction이 제거하는 편향 항과 남기는 항을 구분한다. Bagging identity만으로 원래 target과의 동일성을 주장하지 않는다.
3. 타당한 점추정량을 정한 뒤 그 추정량의 표준오차와 단측 추론을 검토한다. Nested bootstrap을 그 자체로 해결책으로 취급하지 않는다.
4. 이 연결을 predictor weights가 학습되는 ERA/GSCA로 확장하고, 기존 6성분 manuscript simulation과 tuning 절차는 별도로 검증한다.

이번 결과는 fixed-indicator target과 점추정 보정 후보를 평가하는 단계다. GSCA 전체에 대한 finite-$N$ 보정 및 검정의 이론적 타당성은 아직 확립되지 않았다.

## 재현과 결과 파일

- [정의와 API 설명](README_indicator.md)
- [전체 기본 실험 요약 CSV](results_indicator/GSCA_Indicator_Primary_Summary.csv)
- [전체 요약 JSON](results_indicator/GSCA_Indicator_Primary_Results.json)
- [코드·원자료별 결과·참값·seed 재현 패키지](results_indicator/GSCA_Indicator_Primary_Package.zip)
- [기본 실험 실행 기록](https://github.com/PsycheMatrica/GSCA.Reg_Prime/actions/runs/42336898888)
- [검증 및 ACSI 실행 기록](https://github.com/PsycheMatrica/GSCA.Reg_Prime/actions/runs/42337087534)

기본 실험의 네 primary job은 모두 성공했다. 같은 run의 미완료 nested job만 이후 중단했기 때문에 run 전체 상태는 cancelled로 표시된다. 이 보고서는 완료된 primary 결과만 사용한다. 기본 통계 계산 commit은 c0f48f19c7c4f4aab2aa6f5b7d5d9572be453fa2이며, ACSI 검증 commit은 52691d91cc829d7119fc1c56fbac17a93c679e4d이다.
