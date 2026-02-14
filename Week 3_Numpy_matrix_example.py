import numpy as np

##### 행렬과 벡터의 곱셈

# 가상의 환자 데이터 (행: 환자 3명, 열: [소득분위, 나이, BMI])
X = np.array([
    [1, 45, 22],
    [2, 32, 28],
    [3, 55, 30]
])

#print(X.shape)

# 종속 변수 (예: 표본 가중치, 회귀계수)
y = np.array([80, 0.5, 1.2])

#print(y.shape)    # 행이 3개인 1차원 벡터

# 행렬 X와 벡터 y의 요소별 곱
# z=X*y
# print(z)   # y가 가중치일 때 유용

# 행렬 X와 벡터 y의 행렬 곱(앞 행렬의 열 개수와 뒤 행렬의 행 개수가 반드시 같아야)
z=X@y
print(z)      
print("환자별 예측 혈압:", z)

# 전치
V=X.T
print(V)
print(V.shape)

# 브로드캐스팅 덧셈
X2 = np.array([[1, 45, 22]]).T
print(X2.shape)

W = np.array([[10, 20]])
print(W.shape)

sum = X2 + W
print(sum)


#####  벡터의 길이(크기) 계산 --> 환자의 특성간 거리에 활용

# 가상의 환자 데이터 (예: 혈압 변화량 벡터)
v = np.array([3, -4, 0])

# 방법 1: NumPy의 내장 함수 사용 (L2 Norm)
norm_value = np.linalg.norm(v)

# 방법 2: 수식을 직접 구현 (요소의 제곱합의 제곱근)
sum_sq = np.sum(v**2)
manual_value = np.sqrt(sum_sq)

print(f"벡터의 길이(원소 개수): {len(v)}")
print(f"내장 함수 결과: {norm_value}")
print(f"직접 계산 결과: {manual_value}")


#####  벡터 집합과 선형 가중 결합

# 1. 두 가지 건강 지표 벡터 (환자 3명)
v1 = np.array([2.5, 2, 1.5]) # 운동 시간 (시간/주)
v2 = np.array([5, 3, 8]) # 하루 흡연량 (개비)

# 2. 각 지표에 부여할 가중치(선형회귀계수)
w1 = 10  # 운동의 긍정적 연관성 
w2 = -5  # 흡연의 부정적 연관성

# 3. 선형가중결합 계산 (w1*v1 + w2*v2)  - 여러 변수의 영향력을 가중치로 조절하여 하나의 통합된 지표로 표현
health_score = (w1 * v1) + (w2 * v2)

print("운동 지표 (v1):", v1)
print("흡연 지표 (v2):", v2)
print("건강 수준 점수:", health_score)   


##### 코사인 유사도(Cosine Similarity)와 피어슨 상관계수(Pearson Correlation) 

# 다섯 명의 환자에 대한 나이와 BMI
x = np.array([40, 20, 68, 32, 82])
y = np.array([25.3, 24, 22.4, 18, 17.5])

# 각 벡터의 평균 구하기
mean_x = np.mean(x)
mean_y = np.mean(y)
print(mean_x)

#  데이터 센터링: 각 값에서 평균을 빼기
x_centered = x - mean_x
y_centered = y - mean_y
print(x_centered)

#  두 벡터의 내적(Dot Product)
dot_product = np.dot(x_centered, y_centered)

#  분모: 각 벡터의 노름(Norm)의 곱
norm_x = np.linalg.norm(x_centered)
norm_y = np.linalg.norm(y_centered)
magnitude = norm_x * norm_y

# 코사인 유사도 계산
cosine_sim = dot_product / magnitude

# 피어슨 상관계수와 비교 (Numpy 내장 함수)
# [0, 1] 위치의 값이 x와 y 사이의 상관계수
pearson_corr = np.corrcoef(x, y)[0, 1]

print(f"계산된 코사인 유사도: {cosine_sim:.4f}")   # 나이와 BMI간의 상관관계(유사도)
print(f"피어슨 상관계수: {pearson_corr:.4f}")

##### Numpy로 k-means clustering 구현하기

# 1. 중심점 초기화: 데이터 중 임의로 K개의 점을 선택해 각 군집의 중심으로 정하기
# 2. 군집 할당: 모든 데이터 포인트에서 각 중심점까지의 거리를 계산해, 가장 가까운 중심점의 군집으로 할당
# 3. 중심점 위치 조정: 각 군집에 속한 데이터들의 평균값을 구해 새로운 중심점을 정하고 변화가 없을 때까지 이 과정을 반복

# 1. 중심점 초기화
# 8명 환자, 3개 그룹의 [나이, 수축기 혈압] 자료
X = np.array([
    [25, 110], [28, 115], [30, 120],  # 그룹 A 
    [60, 140], [65, 145], [70, 150],  # 그룹 B 
    [45, 130], [50, 135]              # 그룹 C 
])

# 군집 수(K) 설정 및 초기 중심점 무작위 선택
K = 2
indices = np.random.choice(len(X), K, replace=False)  
centroids = X[indices]

print("초기 중심점:\n", centroids)

# 2. 군집 할당

# 각 데이터(X)에서 각 중심점(centroids)까지의 거리 계산
# (8, 1, 2) 형태와 (2, 2) 형태의 계산을 위해 차원을 조절합니다.
distances = np.linalg.norm(X[:, np.newaxis] - centroids, axis=2)
print(distances)   # 각 환자(행)가 0번 중심점과 1번 중심점에서 각각 얼마나 떨어져 있는지, 첫번째 환자는 두번째 중심점과 더 가까워 두번째 군집으로 분류됨

# 가장 거리가 짧은 중심점의 인덱스 선택 (0 또는 1)
labels = np.argmin(distances, axis=1)

print("각 데이터의 소속 군집:", labels)   # 환자 8명의 가장 가까운 중심점 번호

# 3. 중심점 위치 조정
new_centroids = np.array([X[labels == k].mean(axis=0) for k in range(K)])


print("새로운 중심점:\n", new_centroids)
