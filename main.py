import scipy
import numpy as np
import matplotlib.pyplot as plt

with open("Output.txt") as f:
    detector = f.readline()
    while "LC Chromatogram" not in detector:
        detector = f.readline()
    while "R.Time" not in detector:
        detector = f.readline()
    Time = []
    Intensity = []
    while detector != "\n":
        detector = f.readline()
        if detector != "\n":
            t, i = detector.split()
            Time.append(float(t))
            Intensity.append(float(i))

f = scipy.interpolate.interp1d(Time, Intensity, "cubic")
x_range = np.linspace(2.2, 3.2, 10000)
plt.figure(figsize=(10, 4))
plt.title("Analysis Figure")
plt.xlabel("t/min")
plt.ylabel("A/uV")
plt.plot(x_range, f(x_range))


def maximum(max_min, max_max, func, eps=10 ** -5):
    left, right = 0, 0
    while max_max - max_min > eps:
        left = (max_max - max_min) / 3 + max_min
        right = (max_max - max_min) / 3 * 2 + max_min
        if func(left) > func(right):
            max_max = right
        else:
            max_min = left
    return max(left, right)


def inflection(inf_min, inf_max, func, eps=10 ** -5):
    while inf_max - inf_min > eps:
        points = np.linspace(inf_min, inf_max, 5)
        de_diff = np.diff(np.diff(func(points)))
        if de_diff[0] > 0:
            if de_diff[1] > 0:
                if de_diff[2] > 0:
                    inf_min = points[2]
                else:
                    inf_min = points[1]
            else:
                if de_diff[2] > 0:
                    inf_min = points[2]
                    inf_max = points[3]
                else:
                    inf_max = points[3]
        else:
            if de_diff[1] > 0:
                inf_min = points[1]
                inf_max = points[2]
            else:
                inf_max = points[2]
    return (inf_min + inf_max) / 2


max_1 = maximum(2.52, 2.54, f)
max_2 = maximum(2.86, 2.88, f)
plt.vlines((max_1, max_2), 0, (f(max_1), f(max_2)), linestyles="--")
plt.hlines(f(max_1), 2.3, max_1, linestyles=":")
plt.annotate("h", (2.35, f(max_1)), (2.35, f(max_1) / 2), arrowprops={"arrowstyle": "->"},
             horizontalalignment="center", verticalalignment="center")
plt.annotate("h", (2.35, 0), (2.35, f(max_1) / 2), arrowprops={"arrowstyle": "->"},
             horizontalalignment="center", verticalalignment="center")
a_inf_1 = inflection(2.47, 2.53, f)
d_inf_1 = inflection(2.54, 2.6, lambda x: -f(x))
a_inf_2 = inflection(2.79, 2.86, f)
d_inf_2 = inflection(2.87, 2.96, lambda x: -f(x))
plt.vlines((a_inf_1, a_inf_2), 0, (f(a_inf_1), f(a_inf_2)), linestyles="--")
plt.vlines((d_inf_1, d_inf_2), 0, (f(d_inf_1), f(d_inf_2)), linestyles="--")
plt.scatter((a_inf_1, d_inf_1, a_inf_2), (f(a_inf_1), f(d_inf_1), f(a_inf_2)), 5, "C0", "s")
plt.scatter((d_inf_2, max_1, max_2), (f(d_inf_2), f(max_1), f(max_2)), 5, "C0", "s")
plt.vlines((a_inf_1, d_inf_1), (f(a_inf_1), f(d_inf_1)), 200000, linestyles=":")
plt.annotate("", (a_inf_1, 180000), (d_inf_1, 180000), arrowprops={"arrowstyle": "<->"})
plt.text((d_inf_1 + a_inf_1) / 2, 190000, "2σ", horizontalalignment="center", verticalalignment="center")
N_1 = (max_1 / (d_inf_1 - a_inf_1) * 2) ** 2
N_2 = (max_2 / (d_inf_2 - a_inf_2) * 2) ** 2
print("N_1 = " + str(N_1))
print("N_2 = " + str(N_2))
print()

points_005h_1 = scipy.optimize.fsolve(lambda x: f(x) - 0.05 * f(max_1), np.array((2.46, 2.62)))
points_005h_2 = scipy.optimize.fsolve(lambda x: f(x) - 0.05 * f(max_2), np.array((2.79, 2.98)))
plt.hlines(0.05 * f(max_1), points_005h_1[0], points_005h_1[1], linestyles="--")
plt.hlines(0.05 * f(max_2), points_005h_2[0], points_005h_2[1], linestyles="--")
plt.scatter(points_005h_1, [0.05 * f(max_1)] * 2, 5, "C0", "s")
plt.scatter(points_005h_2, [0.05 * f(max_2)] * 2, 5, "C0", "s")
plt.vlines(np.append(points_005h_1, max_1), [0.05 * f(max_1)] * 2 + [210000], 230000, "C0", ":")
plt.hlines(0.05 * f(max_1), 2.4, points_005h_1[0], linestyles=":")
plt.annotate("0.05h", (2.41, 0.05 * f(max_1)), (2.41, 50000), arrowprops={"arrowstyle": "->"},
             horizontalalignment="center", verticalalignment="center")
plt.annotate("", (points_005h_1[0], 220000), (max_1, 220000), arrowprops={"arrowstyle": "<->"})
plt.annotate("", (max_1, 220000), (points_005h_1[1], 220000), arrowprops={"arrowstyle": "<->"})
plt.text((points_005h_1[0] + max_1) / 2, 230000, "A", horizontalalignment="center", verticalalignment="center")
plt.text((max_1 + points_005h_1[1]) / 2, 230000, "B", horizontalalignment="center", verticalalignment="center")
fs_1 = (points_005h_1[1] - points_005h_1[0]) / (max_1 - points_005h_1[0]) / 2
fs_2 = (points_005h_2[1] - points_005h_2[0]) / (max_2 - points_005h_2[0]) / 2
print("fs_1 = " + str(fs_1))
print("fs_2 = " + str(fs_2))
print()

R = 2 * (max_2 - max_1) / (2 * (d_inf_1 - a_inf_1) + 2 * (d_inf_2 - a_inf_2))
print("R = " + str(R))

plt.show()
