# GSCA의 OOB 예측 기여 보정: 구현과 탐색적 시뮬레이션

작성일: 2026-11-06 (미국 동부시간)

## 목적과 현재 판단

관심 모수는 정규화하지 않은 $\Delta\mathrm{EPE}_p(N)$이다. Predictor replacement는 해당 성분을 같은 식의 다른 부모 성분에 대한 선형투영으로 대체한다. 모델은 전체모형을 한 번 적합하며, 제거 후 outcome 식을 다시 적합하지 않는다.

이번 작업은 기존 GSCA 코드의 연구용 확장과 검증용 simulation v3를 추가한 것이다. 기존 공개 함수의 기본 동작은 유지한다. 기존 6성분 ACSI simulation 전체를 재현한 결과가 아니라, 원인을 구분할 수 있는 3성분 모형에서의 단계적 검증이다.

**핵심 구분은 점추정의 편향 보정과 단측 percentile 검정의 보정이다.** 점추정이 개선되었다고 같은 보정량만큼 percentile 분포를 이동하면 검정도 적절해지는 것은 아니다.

## 왜 잔차 분산을 보정하는가

GSCA 구조모형에서 평가하는 내생 성분 $t_b$는 bootstrap마다 추정 가중치와 함께 변한다. 전체 예측을 $f_b$, 제거 후 예측을 $g_{p,b}$, 제거되는 항을 $h_{p,b}=f_b-g_{p,b}$라 하면,

$$
C_{p,i}
=s^2(t_b-f_b)-s^2(t_b-g_{p,b})
=s_f^2-s_g^2-2\widehat{\operatorname{Cov}}(t_b,h_{p,b}).
$$

따라서 예측값의 분산 차이만 사용하면 내생 성분과 제거되는 항의 공분산을 놓친다. 코드에서는 paired residuals를 직접 저장·집계한다.

관측값 $i$가 OOB인 반복 수를 $K_i$라고 할 때,

$$
\widehat{\Delta\mathrm{EPE}}_p^{\mathrm{raw}}
=\frac1N\sum_i\overline{e_{-p,bi}^2-e_{1,bi}^2},
\qquad
\widehat{\Delta\mathrm{EPE}}_p^{\mathrm{VC}}
=\widehat{\Delta\mathrm{EPE}}_p^{\mathrm{raw}}
+\frac1N\sum_i\{s^2(e_{1,bi})-s^2(e_{-p,bi})\}.
$$

여기서 각 관측값은 동일 가중을 받고, 분산의 분모는 $K_i-1$이다. 보정량은 항상 양수일 필요가 없다.

다음은 정확한 계산 항등식이다.

$$
\widehat{\Delta\mathrm{EPE}}_p^{\mathrm{VC}}
=\frac1N\sum_i(\bar e_{-p,i}^2-\bar e_{1,i}^2)
+\frac1N\sum_i\frac{s^2(e_{1,bi})-s^2(e_{-p,bi})}{K_i}.
$$

이 항등식은 bagging 및 유한 bootstrap 반복 수와의 관계를 보여준다. **원래 unbagged GSCA의 유한-$N$ 위험에 대한 불편성을 증명하지는 않는다.** 추정 가중치, penalization 및 비선형 추정의 고차 편향이 남을 수 있으므로 독립적으로 계산한 참값과 비교해야 한다.

## 구현한 단측 방법

| 방법 | percentile에 사용하는 값 |
|---|---|
| 기존 규칙 비교 | 계수 0 대체를 사용한 개별 OOB 반복 손실 차이 |
| Raw percentile | projection 대체를 사용한 개별 OOB 반복 손실 차이 |
| Shifted percentile | raw 반복값에 추정 보정량을 더한 값 |
| Nested raw percentile | 바깥 pairs bootstrap에서 다시 계산한 최종 raw 추정량 |
| Nested corrected percentile | 바깥 pairs bootstrap에서 다시 계산한 최종 VC 추정량 |

모든 predictive bound는 **단측 95% 하한**이다. 양측 predictive CI를 사용하지 않았다. 경험적 $p$-value는 $B^{-1}\sum_b1(T_b\le0)$이며, 하한이 0보다 큰 경우 기각한다. 유한 반복에서의 동률을 일관되게 처리하기 위해 $p<.05$를 사용한다. 0이라는 수치는 nonpositive 반복값을 관찰하지 못했다는 뜻이다.

기존 계수 0 대체와 projection 대체는 서로 다른 모수를 정의한다. 따라서 기존 규칙의 기각률을 projection 경계에서 비교할 수는 있지만, 그것을 기존 방법 자체의 영가설에 대한 Type I error라고 부르면 안 된다.

Nested 절차는 원래 추정량 전체를 재계산한다. 그러나 이것만으로 quadratic degeneracy가 해결되지는 않는다. 이에 관한 일반적 논의는 [Chen and Fang (2019)](https://arxiv.org/abs/1901.04861)를 참고할 수 있다. 이번 GSCA 결과의 타당성은 별도의 simulation 근거로 판단한다.

## 평가 대상과 참값

성분은 $P=3$, 성분별 지표 수는 $J_p=1$ 또는 3이다. 두 외생 성분의 상관은 .30이고, 내생 성분에 대한 nuisance 경로는 .35이다. 관심 경로는 0, predictive boundary, .30으로 나눈다. $J_p=1$은 OLS 연결을 확인하고, $J_p=3$은 가중치와 내생 성분을 실제로 학습한다. 데이터 생성은 성분이 지표의 정확한 선형결합이 되도록 구성했다.

참값은 각 독립 훈련자료로 적합한 affine residual map의 위험을 알려진 모집단 평균·공분산으로 정확히 적분하고, 이를 훈련자료들에 대해 평균하여 계산한다. Test-sample simulation 오차는 없다. 훈련자료 평균의 Monte Carlo 오차는 남으므로 별도로 보고한다.

정규분포·무벌점($\lambda=0$) $J_p=1$ 조건에서는 정확히

$$
\Delta\mathrm{EPE}_p(N)
=.91b^2-\frac{.8775-b^2-.21b}{N-4}
$$

이다. $N=80$의 정확한 양의 경계는 $b=.1103486992$이다. 이 조건의 bias/coverage는 analytic truth로 계산한다. 별도의 훈련자료 simulation은 그 공식을 검증한다. 그 외 GSCA 경계는 Monte Carlo로 보정한 근사 경계이며, 실제 검정 크기 .05를 맞추려고 경로값이나 보정 계수를 조정한 것이 아니다.

Oracle reference는 모집단 공분산을 알고 있는 경우다. Sample reference는 원자료 전체의 지표 공분산을 사용한다. 후자는 내생 성분의 척도 추정도 포함하므로 predictor-only nuisance estimate라고 해석하면 안 된다.

구조모형의 손실은 학습된 내생 성분의 위험이다. Reference variance로 척도를 맞춰도 그 성분의 정체성이 고정되는 것은 아니다. 측정모형의 C-edge 손실은 해당 지표가 성분 구성에 포함될 수 있는 reconstruction 위험이며, 관찰 불가능한 outcome의 forecasting 위험과 구별한다. 이번 simulation의 관심 경로는 B-edge이다.

## Nested percentile의 직접 비교

조건은 $N=80,\ J_p=3,\ \lambda=0$, 정규 자료, oracle reference이며 조건별 원자료 100개를 사용했다. 원자료의 OOB 반복은 199회, outer bootstrap은 99회, 각 outer sample의 inner OOB 반복은 99회다.

| 방법 | 예측 경계 부근 기각률 | 양의 기여 조건 검정력 |
|---|---:|---:|
| Raw percentile | 0.000 | 0.020 |
| Shifted percentile | 0.000 | 0.020 |
| Nested raw percentile | 0.020 | 0.310 |
| Nested corrected percentile | 0.050 | 0.440 |

경계 후보의 독립적인 위험 추정은 $-0.000039$, MCSE는 $0.000372$였다. 정확히 알려진 경계가 아니라, 위험이 거의 0인 것으로 확인된 조건이다. 양의 기여 조건의 위험은 $0.072193$, MCSE는 $0.000316$이었다. 각각 2,000개의 독립 훈련자료로 계산했다.

Nested corrected의 95% 단측 하한 coverage는 경계 부근 .95, 양의 기여 조건 .94였다. 관측 기각률 .05의 Wilson 95% Monte Carlo 구간은 [.0215, .1118], 검정력 .44의 구간은 [.3467, .5377]이다. 이 수치만으로 정확한 검정 크기나 일양적 Type I error 제어를 확립했다고 할 수 없다.

이 비교는 두 역할을 구분해서 살펴본다. Raw에서 nested raw로의 변화는 percentile에 넣는 통계량을 최종 추정량 전체로 바꾸는 효과를 보여주고, nested raw에서 nested corrected로의 변화는 점추정 보정을 추가한 효과를 보여준다. Raw와 nested의 비교에는 재표집 구조와 inner Monte Carlo 반복 수의 차이가 함께 들어가므로 순수한 단일 요인의 효과 분해는 아니다. Nested raw와 corrected는 동일한 outer/inner 표본에 대한 paired 비교다. 관측된 개선은 이번 조건에 대한 탐색적 근거다.

점추정 bias는 경계 조건에서 raw −.011064, VC +.001006이었다. 양의 조건에서는 raw −.007935, VC +.002161이었다. 유한 simulation에서 관찰한 bias이며, 불편성의 정리가 아니다.


## 표본 reference 및 이분산에 대한 추가 비교

$J_p=3,\lambda=0$에서 $N=80,200$, 정규·이분산, 계수 0·경계 후보·양의 기여, oracle·sample reference를 교차한 24개 조건/reference 행을 조사했다. 각 행은 원자료 200개, OOB 199회다. Oracle과 sample reference는 같은 원자료와 bootstrap seed를 사용한 paired 비교이므로 24행이 서로 독립적인 24실험이라는 뜻은 아니다.

**24개 행 모두에서 관측된 절대 bias와 RMSE가 줄었다.** 아래는 실제 자료에 가까운 sample-reference 경계 후보 결과다.

| N | 자료 | Raw bias | VC bias | Raw RMSE | VC RMSE |
|---:|---|---:|---:|---:|---:|
| 80 | 정규 | -0.0125 | 0.0003 | 0.0356 | 0.0319 |
| 80 | 이분산 | -0.0191 | 0.0006 | 0.0564 | 0.0501 |
| 200 | 정규 | -0.0054 | -0.0008 | 0.0113 | 0.0097 |
| 200 | 이분산 | -0.0076 | 0.0005 | 0.0225 | 0.0208 |

이 조건의 경계 참값에는 Monte Carlo 오차가 남아 있다. CSV의 TruthMean, TargetMCSE 및 Classification을 함께 보아야 한다. 경계 후보의 raw/shifted 기각률은 .000–.005였으므로, 점추정 개선이 개별 OOB 반복값의 percentile 검정 보수성을 해결한 것은 아니다. $N=200$, 정규·양의 기여·sample reference에서는 검정력이 .155에서 .190으로 늘었지만, 이 역시 nested 방법의 광범위한 검증을 대신하지 않는다.

초기 12개 행에서는 11개 행에서 절대 bias와 RMSE가 줄었다. 예외는 단일 지표, $\lambda=.1$, 양의 기여 조건이었다. Bias가 −.00126에서 +.00607로 커졌고 RMSE는 .05570에서 .05574로 조금 늘었다. 따라서 모든 조건에서 개선된다고 일반화하지 않는다.

## 실행 범위와 검증

- 최종 실험: 38개 조건/reference 행, 원자료 수준 분석 6,200회. 같은 자료의 oracle/sample 비교를 포함한 분석 횟수다.
- OOB 적합: 3,194,000회. 별도의 원자료 적합, calibration/truth 적합 및 단위 검증은 이 수에 포함하지 않았다.
- 분석된 원자료의 적합 실패 및 OOB 적합 실패: 0회. Nested 두 조건도 각각 100/100개 자료에서 outer inference를 계산했다.
- MATLAB R2026a에서 OLS/ridge 연결, 측정모형만 있는 경우, 성분 순서 변경, 최종 목적함수, 수렴 진단, singular Gram 및 부호 anchor 처리, affine risk 계산, projection 직교성, 척도·부호 변환, paired residual 항등식, nested 실행을 검사했다.
- 수치 overflow와 유한 반복 percentile grid의 경계 처리도 보완하고 별도로 검증했다.
- Gaussian OLS의 analytic truth와 독립 training-risk 계산은 세 효과 조건에서 모두 1 MCSE 이내로 일치했다.

수치가 0회 실패였다는 사실이 모든 모집단에서 안정성이나 위험의 적분 가능성을 증명하지는 않는다. 코드는 실패를 자동으로 다시 뽑아 숨기지 않는다.

## 연구상 판단과 다음 검증

현재 가장 유망한 후속 검증 후보는 **paired residual variance correction + 전체 추정량의 nested one-sided percentile**이다. 단순한 shifted percentile은 점추정 bias 개선을 검정 개선으로 충분히 연결하지 못했다.

다만 현재 nested 성능은 정규, 무벌점, $N=80$, oracle reference, 조건별 100개 자료와 outer/inner 각각 99회에 국한된다. 다음 항목이 확인되기 전에는 검증된 실무 검정으로 제시하지 않는다.

1. 표본 reference, 이분산, 큰 N, 고정 양의 벌점에서 nested 검정의 경계 기각률과 검정력.
2. 원자료 반복 수 및 outer/inner bootstrap 수를 늘렸을 때의 안정성. 특히 inner Monte Carlo noise가 관측 기각률을 높였는지 구분.
3. Monte Carlo 경계를 더 정밀하게 확인하고, 영가설 내부 및 경계 양쪽의 가까운 조건을 조사.
4. 기존 6성분 ACSI simulation과 자료 예시로 확장.
5. $\lambda$를 선택하는 알고리즘까지 목표에 포함하려면 in-bag 내부에서 tuning을 반복하는 구현.
6. 일반 GSCA에서 원래 유한-$N$ target에 대한 보정 전개와 경계 추론의 이론적 정당화.

현재 GitHub 원본은 이미 EPE 차이를 계산하므로, 코드의 EPE 분모를 새로 삭제할 필요는 없었다. 이번 변경의 핵심은 replacement, 평가 척도 명시, 잔차 기반 보정 및 percentile에 넣는 통계량이다.

## 재현 및 파일

```matlab
addpath('research');
test_gsca_fit;
test_gsca_risk_maps;
test_gsca_predictive;

run_research_experiment('primary');
run_research_experiment('sensitivity');
run_research_experiment('nested_boundary');
run_research_experiment('nested_positive');
```

- [연구 코드 및 방법 설명](README.md)
- [전체 요약 CSV](results/GSCA_Correction_Summary.csv)
- [전체 결과와 참값 JSON](results/GSCA_Correction_Results.json)
- [코드·CSV·MAT 재현 패키지](results/GSCA_Correction_Research_Package.zip)
- [통계 실행 기록](https://github.com/PsycheMatrica/GSCA.Reg_Prime/actions/runs/42323844916)
- [최종 코드 검증 기록](https://github.com/PsycheMatrica/GSCA.Reg_Prime/actions/runs/42324107037)

통계 실행 commit은 a2b16cb509e51bb447b975c8ecfec64f97316a12, 추가 수치 방어와 percentile grid 처리를 검증한 코드 commit은 b9546acbeea42a0fcc7929abdd5ebe0ecf8ea45e이다. 후자의 수정은 이번 유한 수치와 $\alpha=.05$, 99/199회 설정의 계산값을 바꾸지 않는다. 패키지는 두 버전을 manifest에 명시하며 원자료별 결과, 설정 및 재현 코드를 포함한다.
