"""Draw the implemented one-joint Hill system, not an anatomical reconstruction.

Geometry and displayed numerical values come from the portable calibration.
Muscle paths/limb outlines are illustrative; the solver uses constant effective
moment arms and linear MT length maps. No independent BSE mechanism or assist
torque is present in this demo. PNG previews and editable SVG/PDF are exported.
"""
from pathlib import Path
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/hill_compliant_evaluation"
BLUE = "#1764ad"
ORANGE = "#d65b22"
PURPLE = "#7952a8"
DARK = "#24354a"
GREY = "#79879a"
TEAL = "#18847a"

LABELS = {
    "zh": {
        "title": "当前 BSE 项目的等效肌肉—骨骼模型",
        "subtitle": "固定骨盆 · 单转动自由度 · 两个 Hill 肌肉单元 · 弹性肌腱",
        "a": "A   身体、负载与相反方向的肌肉力矩",
        "b": "B   每个肌肉单元的内部结构",
        "c": "C   同一组 Hill 力驱动角度并输出力矩",
        "trunk": "等效刚性躯干",
        "com": "质心",
        "joint": "O：等效 L5/S1 关节",
        "back": "背部等效伸肌",
        "front": "前部等效屈肌",
        "load": "负载",
        "horizontal": "水平参考线（与地面平行）",
        "ground": "固定骨盆 / 下肢",
        "soft": "关节软限位",
        "damping": "关节速度阻尼",
        "ce": "CE  主动收缩",
        "pe": "PE  被动纤维弹性",
        "fiber_damping": "纤维阻尼",
        "tendon": "SE  弹性肌腱",
        "neural": "激活",
        "two": "背部、前部各一套相同结构；各自参数可标定",
        "parameters": "全九组标定的共享身体参数",
        "foot": "示意姿态：theta=0.95 rad。走线/肢体轮廓为示意，力箭头未按大小比例绘制。",
        "assist": "当前 demo 未加入独立外骨骼助力输入；软限位单独列出，不并入 Hill torque。",
        "net": "净 Hill 力矩",
        "soft_note": "软限位由预测角度越界产生，不追踪 reference。\n图示姿态位于区间内，因此此时软限位力矩为零。",
        "flow_title": "同一 Hill 系统如何产生 torque 与 theta",
        "experiment": "实验来源输入\n背部 + 采集的前部 activation",
        "map": "MAP 来源输入\n背部 + 辅助生成的前部 activation",
        "choose": "每次运行选择一种输入来源；模型形式一致",
        "activation": "激活动力学 + 估计募集\n两种来源可分别标定参数",
        "geometry": "肌肉几何\n由预测角度更新肌肉总长",
        "pair": "两个 Hill + 弹性肌腱\n联立纤维 / 肌腱力平衡",
        "joint_block": "躯干运动平衡\nHill + 重力 + 阻尼 + 软限位",
        "output": "输出\n净 Hill torque 和 theta",
        "feedback": "预测角度更新肌肉几何\n同一步内联立求解",
        "known": "已知负载 m(t)\n共享的 M, L1, L2, phi\n关节阻尼与软限位",
        "reference": "reference 用于训练与评价；只在开始时提供 theta0、omega0。\n预测过程中不读取后续 reference，不按它反算 torque。",
    },
    "en": {
        "title": "Equivalent musculoskeletal model used in the current BSE demo",
        "subtitle": "Fixed pelvis | One rotational DOF | Two opposing Hill units | Compliant tendons",
        "a": "A   Body, load and opposing muscle torques",
        "b": "B   Internal structure of each muscle unit",
        "c": "C   The same Hill forces drive motion and produce torque",
        "trunk": "Equivalent rigid trunk",
        "com": "COM",
        "joint": "O: equivalent L5/S1 joint",
        "back": "Lumped back extensor",
        "front": "Lumped front flexor",
        "load": "Load",
        "horizontal": "Horizontal reference (parallel to ground)",
        "ground": "Fixed pelvis / lower limbs",
        "soft": "Joint soft limits",
        "damping": "Joint viscous damping",
        "ce": "CE: active contraction",
        "pe": "PE: passive fiber elasticity",
        "fiber_damping": "Fiber damping",
        "tendon": "SE: compliant tendon",
        "neural": "Activation",
        "two": "The back and front use this structure, with separately calibrated muscle parameters.",
        "parameters": "Shared body geometry from the all-trial calibration",
        "foot": "Illustrative pose: theta=0.95 rad. Paths/outlines are schematic; force arrows are not magnitude-scaled.",
        "assist": "No independent exoskeleton assistance input is applied; soft-limit torque is reported separately.",
        "net": "Net Hill torque",
        "soft_note": "Soft limits use predicted angle only; no reference tracking.\nThe illustrated pose is within the interval: soft-limit torque is zero.",
        "flow_title": "Torque and angle emerge from the same coupled Hill system",
        "experiment": "Experimental inputs\nBack + recorded front activation",
        "map": "MAP inputs\nBack + auxiliary generated front activation",
        "choose": "Choose one source per rollout; identical model equations",
        "activation": "Activation dynamics + estimated recruitment\nSource-specific calibration is allowed",
        "geometry": "Muscle geometry\nMT lengths from the predicted angle",
        "pair": "Two Hill units + compliant tendons\nCoupled fiber / tendon force equilibrium",
        "joint_block": "Trunk motion balance\nHill + gravity + damping + soft limits",
        "output": "Outputs\nNet Hill torque and angle",
        "feedback": "Predicted angle updates geometry\nCoupled within the same implicit step",
        "known": "Known load m(t)\nShared M, L1, L2, phi\nJoint damping and soft limits",
        "reference": "References are used for calibration/evaluation and the initial theta0, omega0 only.\nFuture reference samples do not enter the rollout or generate predicted torque.",
    },
}


def configure(language):
    font = "Microsoft YaHei" if language == "zh" else "Arial"
    try:
        font_manager.findfont(font, fallback_to_default=False)
    except ValueError:
        font = "DejaVu Sans"
    plt.rcParams.update({
        "font.family": font, "font.size": 14, "text.color": DARK,
        "axes.unicode_minus": False, "mathtext.fontset": "dejavusans",
        "svg.fonttype": "none", "pdf.fonttype": 42,
    })


def arrow(ax, start, end, color=DARK, style="-|>", lw=1.8, scale=14, **kwargs):
    patch = FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=scale,
                            linewidth=lw, color=color, **kwargs)
    ax.add_patch(patch)
    return patch


def line(ax, points, **kwargs):
    xy = np.asarray(points)
    ax.plot(xy[:, 0], xy[:, 1], **kwargs)


def spring(ax, start, end, color, amplitude=.11, turns=6, lw=2):
    start, end = np.asarray(start), np.asarray(end)
    direction = end-start
    length = np.linalg.norm(direction)
    normal = np.array([-direction[1], direction[0]])/length
    t = np.linspace(0, 1, 2*turns+3)
    offsets = np.zeros(len(t))
    offsets[1:-1] = amplitude*np.where(np.arange(1, len(t)-1) % 2, 1, -1)
    points = start + t[:, None]*direction + offsets[:, None]*normal
    line(ax, points, color=color, lw=lw)


def box(ax, xy, width, height, text="", color=DARK, fill="#f4f7fb", fontsize=14):
    patch = FancyBboxPatch(xy, width, height, boxstyle="round,pad=0.06,rounding_size=0.10",
                           facecolor=fill, edgecolor=color, lw=1.6)
    ax.add_patch(patch)
    if text:
        ax.text(xy[0]+width/2, xy[1]+height/2, text, ha="center", va="center",
                fontsize=fontsize, color=color, linespacing=1.6)
    return patch


def callout(ax, text, position, target, color=DARK, fontsize=14, ha="left"):
    ax.annotate(text, target, xytext=position, ha=ha, va="center", fontsize=fontsize,
                color=color, linespacing=1.5,
                arrowprops={"arrowstyle": "-", "color": color, "lw": 1.2})


def angle_arc(ax, radius, start, end, label, color=GREY, label_radius=None):
    angles = np.linspace(start, end, 80)
    line(ax, np.c_[radius*np.cos(angles), radius*np.sin(angles)], color=color, lw=1.5)
    middle = (start+end)/2
    label_radius = label_radius if label_radius is not None else radius+.04
    ax.text(label_radius*np.cos(middle), label_radius*np.sin(middle), label,
            fontsize=19, color=color, ha="center", va="center")


def dimension(ax, start, end, offset, label, color=GREY):
    start, end, offset = np.asarray(start), np.asarray(end), np.asarray(offset)
    line(ax, [start, start+offset], color=color, lw=.9, ls=":")
    line(ax, [end, end+offset], color=color, lw=.9, ls=":")
    arrow(ax, start+offset, end+offset, color, "<->", lw=1.3, scale=11)
    delta = end-start
    middle = (start+end)/2+offset*1.35
    rotation = np.rad2deg(np.arctan2(delta[1], delta[0]))
    ax.text(*middle, label, color=color, ha="center", va="center", fontsize=17,
            rotation=rotation, bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.8})


def save(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".svg", ".pdf"):
        fig.savefig(OUT / (stem+suffix), dpi=200, facecolor="white")
    plt.close(fig)


def anatomy(parameters, language):
    p, t = parameters, LABELS[language]
    fig = plt.figure(figsize=(19, 11), facecolor="white")
    fig.text(.035, .95, t["title"], fontsize=27, weight="bold")
    fig.text(.035, .91, t["subtitle"], fontsize=17, color=GREY)
    fig.text(.035, .864, t["a"], fontsize=18, weight="bold")
    ax = fig.add_axes([.035, .21, .47, .63])
    ax.set_aspect("equal")
    ax.set_xlim(-.39, 1.0)
    ax.set_ylim(-.50, .73)
    ax.axis("off")
    theta, mass, l1, l2, phi = .95, p[18], p[19], p[12], p[13]
    e = np.array([np.cos(theta), np.sin(theta)])
    n = np.array([-e[1], e[0]])
    el = np.array([np.cos(theta-phi), np.sin(theta-phi)])
    nl = np.array([-el[1], el[0]])
    origin, com, top, load = np.zeros(2), l1*e, 2*l1*e, l2*el
    # Lower-limb contours are a fixed support, not additional model DOFs.
    line(ax, [[-.09, -.07], [.04, -.14], [.14, -.08]], lw=9, color="#d1d8e2")
    line(ax, [[-.07, -.09], [-.09, -.26], [-.16, -.425]], color=GREY, lw=5)
    line(ax, [[.06, -.12], [.11, -.28], [.16, -.425]], color=GREY, lw=5)
    line(ax, [[-.17, -.425], [-.09, -.425]], color=GREY, lw=5)
    line(ax, [[.15, -.425], [.23, -.425]], color=GREY, lw=5)
    line(ax, [[-.34, -.44], [.94, -.44]], color="#d9e0e8", lw=2)
    ax.text(.23, -.48, t["ground"], fontsize=13, color=GREY, va="top")
    line(ax, [origin, top], color="#cbd5e1", lw=24, solid_capstyle="round")
    line(ax, [origin, top], color=DARK, lw=5, solid_capstyle="round")
    ax.add_patch(Circle(top+.075*e, .047, fc="#edf1f6", ec=GREY, lw=2))
    # The drawn arm illustrates the load location; it has no separate state.
    elbow = top+.52*(load-top)+np.array([0, .05])
    line(ax, [top, elbow, load], color="#b5c0ce", lw=5, solid_capstyle="round")
    line(ax, [origin, load], color=TEAL, lw=1.6, ls="--")
    line(ax, [origin, [.92, 0]], color="#c3cbd7", lw=1.2, ls="--")
    horizontal_position = (.73, -.09) if phi > .7 else (.48, -.038)
    ax.text(*horizontal_position, t["horizontal"], fontsize=11, color=GREY, ha="center")

    back_start = np.array([-.135, -.065])
    back_end = .85*top+.065*n
    front_start = np.array([.145, -.065])
    front_end = .85*top-.065*n
    for start, end, color in ((back_start, back_end, BLUE), (front_start, front_end, ORANGE)):
        delta = end-start
        normal = np.array([-delta[1], delta[0]])/np.linalg.norm(delta)
        line(ax, [start, end], color=color, lw=2.5)
        center = start+.57*delta
        ax.add_patch(Ellipse(center, np.linalg.norm(delta)*.45, .045,
                            angle=np.rad2deg(np.arctan2(delta[1], delta[0])),
                            fc=color, ec="white", lw=.8, alpha=.88))
        spring(ax, start+.80*delta, start+.97*delta, color, amplitude=.012, turns=4, lw=1.8)
        ax.add_patch(Circle(start, .009, fc=color, ec="none"))
    callout(ax, t["back"]+"\n"+r"$+R_bF_{T,b}$",
            (-.365, .46), back_start+.61*(back_end-back_start), BLUE, 15)
    callout(ax, t["front"]+"\n"+r"$-R_fF_{T,f}$",
            (.55, -.22) if phi > .7 else (.47, -.145), front_start+.54*(front_end-front_start), ORANGE, 15)
    callout(ax, t["trunk"], (.52, .635), .91*top, DARK, 15)
    ax.add_patch(Circle(com, .017, fc="#efbf46", ec=DARK, lw=1.3, zorder=5))
    ax.text(com[0]+.035, com[1]+.035, "M "+t["com"], fontsize=14)
    arrow(ax, com, com+[0, -.17], DARK, lw=2)
    ax.text(com[0]+.044, com[1]-.135, r"$Mg$", ha="left", fontsize=17)
    ax.add_patch(Rectangle(load+[-.048, -.040], .096, .080, fc="#edf1f6", ec=DARK, lw=1.8))
    ax.text(load[0]+.065, load[1]+.03, t["load"]+"\n"+r"$m(t)$: 5 / 10 / 15 kg",
            fontsize=13, va="center", linespacing=1.6)
    arrow(ax, load+[0, -.04], load+[0, -.215], DARK, lw=2)
    ax.text(load[0]+.025, load[1]-.16, r"$mg$", fontsize=18)
    dimension(ax, origin, com, .10*n, r"$L_1$")
    dimension(ax, origin, load, -.065*nl, r"$L_2$", TEAL)
    angle_arc(ax, .145, 0, theta, r"$\theta$", DARK, .205)
    angle_arc(ax, .445, theta-phi, theta, r"$\phi$", GREY, .49)
    ax.add_patch(Circle(origin, .027, fc="white", ec=DARK, lw=2.3, zorder=6))
    ax.text(-.005, .037, "O", fontsize=17, ha="right")
    callout(ax, t["joint"], (-.36, -.29), origin, DARK, 13)
    # Both purple and gray effects act at the joint. The soft spring is virtual.
    callout(ax, t["soft"]+r" $\tau_{soft}(\theta)$"+"\n"+t["damping"]+r" $-B\omega$",
            (.39, -.325), [.025, -.025], PURPLE, 13)
    ax.add_patch(Circle(origin, .062, fc="none", ec=PURPLE, lw=1.2, ls=":"))
    fig.text(.05, .166, t["parameters"], fontsize=14, color=GREY)
    fig.text(.05, .128, f"M = {mass:.2f} kg     L1 = {l1:.4f} m", fontsize=17)
    fig.text(.05, .093, f"L2 = {l2:.4f} m     phi = {phi:.4f} rad", fontsize=17)

    hill = fig.add_axes([.535, .48, .44, .365])
    hill.set_xlim(0, 10)
    hill.set_ylim(0, 6.6)
    hill.axis("off")
    hill.text(0, 6.35, t["b"], fontsize=18, weight="bold")
    hill.text(0, 5.8, t["two"], fontsize=12, color=GREY)
    left, right = 1.0, 5.2
    line(hill, [[left, 1.7], [left, 4.7]], color=DARK, lw=1.9)
    line(hill, [[right, 1.7], [right, 4.7]], color=DARK, lw=1.9)
    line(hill, [[left, 4.7], [1.8, 4.7]], color=BLUE, lw=2)
    line(hill, [[4.4, 4.7], [right, 4.7]], color=BLUE, lw=2)
    box(hill, (1.8, 4.28), 2.6, .84, t["ce"], BLUE, "#eaf3fc", 14)
    hill.text(.55, 5.25, r"$a_i$", color=BLUE, fontsize=18)
    arrow(hill, [.65, 5.2], [1.78, 4.85], BLUE, scale=10, lw=1.3)
    line(hill, [[left, 3.25], [1.7, 3.25]], color=ORANGE, lw=2)
    spring(hill, [1.7, 3.25], [4.4, 3.25], ORANGE, amplitude=.20)
    line(hill, [[4.4, 3.25], [right, 3.25]], color=ORANGE, lw=2)
    hill.text(3.05, 2.71, t["pe"], fontsize=13, color=ORANGE, ha="center")
    line(hill, [[left, 1.7], [2.7, 1.7]], color=PURPLE, lw=2)
    line(hill, [[2.7, 1.7], [3.4, 1.7]], color=PURPLE, lw=2)
    line(hill, [[3.4, 1.4], [3.4, 2.0]], color=PURPLE, lw=2)
    line(hill, [[2.95, 2.08], [4.05, 2.08], [4.05, 1.32], [2.95, 1.32]], color=PURPLE, lw=2)
    line(hill, [[4.05, 1.7], [right, 1.7]], color=PURPLE, lw=2)
    hill.text(3.2, 1.0, t["fiber_damping"]+r"  $\beta\nu_i$", ha="center", fontsize=13, color=PURPLE)
    line(hill, [[right, 3.55], [6.05, 3.55]], color=DARK, lw=2)
    hill.text(5.60, 3.00, r"$\cos\alpha_i$", fontsize=15, ha="center")
    spring(hill, [6.05, 3.55], [8.50, 3.55], TEAL, amplitude=.20)
    hill.text(7.30, 4.20, t["tendon"], color=TEAL, fontsize=15, ha="center")
    arrow(hill, [8.50, 3.55], [9.75, 3.55], TEAL, lw=2)
    hill.text(9.25, 4.10, r"$F_{T,i}$", fontsize=20, color=TEAL, ha="center")
    hill.text(5, .35, r"$L_{MT,i}=L_{M,i}\cos\alpha_i+L_{T,i}$",
              fontsize=19, ha="center")
    hill.text(5, -.15, r"$F_{T,i}=F_{0,i}\,[a_i f_l f_v+f_{PE}+\beta\nu_i]\cos\alpha_i$",
              fontsize=17, ha="center", clip_on=False)

    balance = fig.add_axes([.535, .115, .44, .30])
    balance.set_xlim(0, 10)
    balance.set_ylim(0, 6)
    balance.axis("off")
    balance.text(0, 5.80, t["c"], fontsize=17, weight="bold")
    box(balance, (.10, 4.43), 9.45, .85,
        t["net"]+r" : $\tau_{Hill}=R_bF_{T,b}-R_fF_{T,f}$", BLUE, "#edf5fc", 18)
    balance.text(.35, 3.65, r"$\tau_g=-MgL_1\cos\theta-mgL_2\cos(\theta-\phi)$", fontsize=18)
    balance.text(.35, 2.66, r"$I(m)\dot{\omega}=\tau_{Hill}+\tau_g-B\omega+\tau_{soft}(\theta)$", fontsize=18)
    balance.text(.35, 1.67, r"$\dot{\theta}=\omega,\qquad I(m)=mL_2^2+\frac{4}{3}ML_1^2$", fontsize=18)
    balance.text(.35, .30, t["soft_note"], color=PURPLE, fontsize=12, linespacing=1.5)
    fig.text(.035, .05, t["foot"], fontsize=13, color=GREY)
    fig.text(.035, .024, t["assist"], fontsize=13, color=PURPLE)
    save(fig, "musculoskeletal_model_"+language)


def flow(language):
    t = LABELS[language]
    fig = plt.figure(figsize=(18, 9), facecolor="white")
    ax = fig.add_axes([.04, .09, .93, .77])
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 10)
    ax.axis("off")
    fig.text(.045, .94, t["flow_title"], fontsize=25, weight="bold")
    fig.text(.045, .89, t["choose"], fontsize=14, color=GREY)
    box(ax, (2.65, 8.55), 4.60, 1.10, t["experiment"], BLUE, "#edf5fc", 14)
    box(ax, (8.50, 8.55), 4.60, 1.10, t["map"], ORANGE, "#fff3ec", 14)
    arrow(ax, [4.95, 8.5], [7.7, 7.45], BLUE)
    arrow(ax, [10.80, 8.5], [8.3, 7.45], ORANGE)
    box(ax, (4.65, 6.55), 6.45, 1.08, t["activation"], DARK, "#f4f7fb", 14)
    arrow(ax, [7.85, 6.49], [7.85, 5.18], DARK)
    box(ax, (.15, 3.30), 4.15, 1.80, t["geometry"]+"\n"+r"$L_{MT,i}(\theta)$", TEAL, "#eef9f7", 14)
    box(ax, (5.20, 3.30), 5.10, 1.80, t["pair"]+"\n"+r"$F_{T,b},\ F_{T,f}$", BLUE, "#edf5fc", 14)
    box(ax, (11.30, 3.30), 4.35, 1.80, t["joint_block"]+"\n"+r"$\theta,\ \omega$", PURPLE, "#f6f0fc", 14)
    arrow(ax, [4.36, 4.15], [5.14, 4.15], TEAL)
    arrow(ax, [10.37, 4.15], [11.24, 4.15], BLUE)
    ax.text(10.83, 4.70, r"$\tau_{Hill}$", fontsize=15, ha="center", color=BLUE)
    box(ax, (11.30, 6.35), 4.35, 1.36, t["known"], GREY, "#f7f8fa", 13)
    arrow(ax, [13.48, 6.29], [13.48, 5.16], GREY)
    # Geometry closure is model feedback, not measured-reference feedback.
    line(ax, [[13.48, 3.24], [13.48, 2.2], [2.2, 2.2], [2.2, 3.24]], color=TEAL, lw=1.9)
    arrow(ax, [2.2, 2.9], [2.2, 3.24], TEAL)
    ax.text(2.30, 1.45, t["feedback"], fontsize=12, color=TEAL, ha="center", linespacing=1.5)
    box(ax, (5.7, .1), 4.7, 1.27, t["output"], DARK, "#f4f7fb", 15)
    arrow(ax, [7.75, 3.23], [7.75, 1.42], BLUE)
    arrow(ax, [13.48, 2.2], [10.45, .85], PURPLE)
    fig.text(.05, .05, t["reference"], fontsize=13, color=GREY, linespacing=1.6)
    save(fig, "model_signal_flow_"+language)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--language", choices=["zh", "en", "both"], default="both")
    args = parser.parse_args()
    data = loadmat(ROOT / "data/processed/hill_compliant/calibrated_models.mat", simplify_cells=True)
    parameters = np.asarray(data["models"][0]["parameters"])
    common = np.array([12, 13, 18, 19, 28, 29, 30, 31])
    assert np.array_equal(parameters[common], np.asarray(data["models"][1]["parameters"])[common])
    for language in (["zh", "en"] if args.language == "both" else [args.language]):
        configure(language)
        anatomy(parameters, language)
        flow(language)
    description = {
        "model": "one-DOF equivalent trunk; two opposing compliant Hill units",
        "musclePaths": "illustrative; solver uses constant moment arms",
        "externalAssistance": "no independent BSE assistance input in this demo",
        "thetaIllustrationRad": .95,
        "sharedAllTrialGeometry": dict(zip(["M_kg", "L1_m", "L2_m", "phi_rad"],
                                          parameters[[18, 19, 12, 13]].tolist())),
        "softLimits": dict(zip(["k_Nm_per_rad", "cap_Nm", "lower_rad", "upper_rad"],
                               parameters[28:32].tolist())),
    }
    (OUT / "model_diagram_metadata.json").write_text(json.dumps(description, indent=2), encoding="utf-8")
    print("Diagrams saved:", OUT)


if __name__ == "__main__":
    main()
