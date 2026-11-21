"""Draw the selected neural/Hill model using the portable parameter cache.

Reuse the existing mechanical illustration; explicitly add input-side neural
feedback. The intent estimate is an excitation input, not the output angle.
"""
import json

import numpy as np
from scipy.io import loadmat

import draw_hill_model as drawing


def flow(language):
    zh = language == "zh"
    title = "当前模型：EMG 神经激励调节 + 同一 Hill 肌肉—躯干系统" if zh else "Selected model: EMG-informed excitation control + coupled Hill mechanics"
    subtitle = "每种来源一组固定参数验证九工况；参考曲线不进入运动中的求解器" if zh else "One fixed vector per source over nine cases; no motion reference enters the solver"
    labels = ([
        "实验来源输入\n背部及前部 activation",
        "MAP 来源输入\n背部 + 辅助前部生成器",
        "已知负载 m(t)\n重力、阻尼、软限位",
        "EMG 意图（补偿负载持续激活）+ 预测姿态反馈\n经同一静态 Hill 方程反求肌肉激励；意图角不是输出角",
        "肌肉几何\n由实际预测角度更新总长",
        "两个 Hill + 弹性肌腱\n纤维 / 肌腱力平衡",
        "躯干运动平衡\n由同一组肌腱力推进状态",
        "预测 theta / omega\n调节肌肉输入",
        "输出：theta、背部力矩、腹部力矩\n净 Hill torque = 背部 − 腹部",
        "角度改变肌肉长度；肌肉力又驱动角度",
        "运动意图由 EMG 历史及已知负载推断。等待通过肌肉激励控制；没有独立外部 PD 力矩。\nreference 保持原 Scale，仅用于离线标定/评价和初始站立姿态；本图不是独立生理验证。",
    ] if zh else [
        "Experimental-source input\nBack + front activation",
        "MAP-source input\nBack + auxiliary front generator",
        "Known load m(t)\nGravity, damping, soft limits",
        "Load-aware EMG intent + predicted-posture feedback\nStatic Hill inversion requests excitation; intent is NOT the output angle",
        "Muscle geometry\nMT lengths from predicted angle",
        "Two Hill + compliant tendons\nFiber / tendon equilibrium",
        "Trunk motion balance\nStates driven by the SAME tendon forces",
        "Predicted angle / speed\nmodify excitation",
        "Outputs: angle, back torque, abdominal torque\nNet Hill torque = back - abdomen",
        "Angle changes lengths; muscle forces drive angle",
        "Intent uses EMG history and known load. Standing control acts through muscle excitation; no external PD torque.\nOriginal Scale is used for offline calibration/evaluation and initial standing pose; no physiological validation is implied.",
    ])
    fig = drawing.plt.figure(figsize=(18, 10), facecolor="white")
    ax = fig.add_axes([.035, .13, .93, .74])
    ax.set_xlim(0, 16.7)
    ax.set_ylim(0, 10)
    ax.axis("off")
    fig.text(.045, .95, title, fontsize=23, weight="bold")
    fig.text(.045, .90, subtitle, fontsize=14, color=drawing.GREY)
    drawing.box(ax, (.2, 8.4), 4.8, 1.2, labels[0], drawing.BLUE, "#edf5fc", 14)
    drawing.box(ax, (5.7, 8.4), 4.8, 1.2, labels[1], drawing.ORANGE, "#fff3ec", 14)
    drawing.arrow(ax, [2.6, 8.32], [3.0, 7.56], drawing.BLUE)
    drawing.arrow(ax, [8.1, 8.32], [7.7, 7.56], drawing.ORANGE)
    drawing.box(ax, (.2, 5.8), 10.3, 1.7, labels[3], drawing.DARK, "#f4f7fb", 14)
    drawing.arrow(ax, [8.1, 5.74], [8.1, 4.48], drawing.DARK)
    ax.text(8.3, 5.0, r"$u_b, u_f$", fontsize=16, color=drawing.DARK)
    drawing.box(ax, (.2, 2.6), 4.8, 1.8, labels[4], drawing.TEAL, "#eef9f7", 14)
    drawing.box(ax, (5.7, 2.6), 4.8, 1.8, labels[5], drawing.BLUE, "#edf5fc", 14)
    drawing.box(ax, (11.2, 2.6), 4.5, 1.8, labels[6], drawing.PURPLE, "#f6f0fc", 13)
    drawing.arrow(ax, [5.08, 3.5], [5.62, 3.5], drawing.TEAL)
    drawing.arrow(ax, [10.58, 3.5], [11.12, 3.5], drawing.BLUE)
    ax.text(10.85, 4.02, r"$\tau_{Hill}$", fontsize=14, color=drawing.BLUE, ha="center")
    drawing.box(ax, (11.2, 7.9), 4.5, 1.45, labels[2], drawing.GREY, "#f7f8fa", 13)
    drawing.arrow(ax, [13.4, 7.82], [13.4, 4.48], drawing.GREY)
    drawing.line(ax, [[11.12, 8.35], [10.8, 8.35], [10.8, 7.1], [10.58, 7.1]], color=drawing.GREY, lw=1.7)
    drawing.arrow(ax, [10.8, 7.1], [10.58, 7.1], drawing.GREY)
    ax.text(11.3, 7.5, "负载持续激活补偿" if zh else "Tonic load compensation", fontsize=11, color=drawing.GREY)
    # Two feedback paths have different purposes: excitation regulation and
    # mechanical geometry closure. Neither reads measured reference states.
    drawing.line(ax, [[15.78, 3.4], [16.25, 3.4], [16.25, 6.45], [10.58, 6.45]], color=drawing.ORANGE, lw=1.8)
    drawing.arrow(ax, [11.05, 6.45], [10.58, 6.45], drawing.ORANGE)
    ax.text(11.25, 6.75, labels[7], fontsize=12, color=drawing.ORANGE)
    drawing.line(ax, [[13.4, 2.52], [13.4, 1.82], [2.6, 1.82], [2.6, 2.52]], color=drawing.TEAL, lw=1.8)
    drawing.arrow(ax, [2.6, 2.2], [2.6, 2.52], drawing.TEAL)
    ax.text(.3, 1.12, labels[9], fontsize=11, color=drawing.TEAL)
    drawing.box(ax, (5.7, .1), 5.4, 1.15, labels[8], drawing.DARK, "#f4f7fb", 13)
    drawing.arrow(ax, [8.1, 2.52], [8.1, 1.33], drawing.BLUE)
    drawing.arrow(ax, [13.4, 1.82], [11.18, .67], drawing.PURPLE)
    fig.text(.045, .05, labels[10], fontsize=12, color=drawing.GREY, linespacing=1.6)
    drawing.save(fig, "model_signal_flow_" + language)


def main():
    model_file = drawing.ROOT / "data/processed/hill_compliant/neural_models.mat"
    package = loadmat(model_file, simplify_cells=True)
    parameters = np.asarray(package["models"][0]["parameters"])
    drawing.OUT = drawing.ROOT / "docs/hill_neural_evaluation"
    for language in ("zh", "en"):
        drawing.configure(language)
        drawing.anatomy(parameters, language)
        flow(language)
    description = {"modelFile": str(model_file.relative_to(drawing.ROOT)),
                   "model": "one-DOF equivalent trunk, two compliant Hill units, input-side neural control",
                   "intentIsOutputAngle": False, "externalPostureTorqueNm": 0,
                   "loadAwareIntent": True,
                   "musclePaths": "Illustrative; constant moment arms used by the solver",
                   "bodyParameters": dict(zip(["M_kg", "L1_m", "L2_m", "phi_rad"], parameters[[18, 19, 12, 13]].tolist())),
                   "softLimits": dict(zip(["stiffness_NmPerRad", "cap_Nm", "lower_rad", "upper_rad"], parameters[28:32].tolist()))}
    (drawing.OUT / "model_diagram_metadata.json").write_text(json.dumps(description, indent=2), encoding="utf-8")
    print("Current model diagrams:", drawing.OUT)


if __name__ == "__main__":
    main()
