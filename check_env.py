"""LS-DYNA + PyDYNA 最小环境验证脚本

用法：
    # 必须使用本项目 venv 中的 Python
    D:\\cxh-daima\\LS-DYNA-ablation\\.venv\\Scripts\\python.exe check_env.py

检查项：
    1) 当前 Python 解释器路径
    2) ansys.dyna.core 是否可导入
    3) PyDYNA 版本（importlib.metadata + 模块自带）
    4) .cursor/rules/pydyna.mdc 是否存在（相对当前工作目录）
    5) shutil.which("lsdyna") 能否找到 LS-DYNA 求解器

附带（不算硬性要求）：
    6) test.k 与 run_test.ps1 是否就位
"""

from __future__ import annotations

import os
import shutil
import sys
from importlib.metadata import PackageNotFoundError, version as pkg_version
from pathlib import Path


# ---- Windows 控制台 UTF-8 设置（避免中文乱码）------------------------------
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass
    try:
        os.system("chcp 65001 > nul")
    except Exception:
        pass


SEP = "=" * 64
RESULTS: list[tuple[str, bool, str]] = []
WARNINGS: list[tuple[str, str]] = []


def record(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))
    flag = "[ OK ]" if ok else "[FAIL]"
    print(f"  {flag}  {name}")
    if detail:
        for line in detail.splitlines():
            print(f"          {line}")


def warn(name: str, detail: str = "") -> None:
    WARNINGS.append((name, detail))
    print(f"  [WARN]  {name}")
    if detail:
        for line in detail.splitlines():
            print(f"          {line}")


def section(title: str) -> None:
    print(f"\n{SEP}\n  {title}\n{SEP}")


# ---------------------------------------------------------------------------
def check_python() -> None:
    section("[1] Python 解释器")
    exe = Path(sys.executable)
    project_root = Path(__file__).resolve().parent
    expected = project_root / ".venv" / "Scripts" / "python.exe"
    is_venv = expected.is_file() and exe.resolve() == expected.resolve()
    ok = is_venv or not expected.exists()
    record(
        "当前 Python 解释器路径",
        ok,
        detail=(
            f"sys.executable = {exe}\n"
            f"项目 venv 路径 = {expected}\n"
            f"sys.version    = {sys.version.split()[0]}"
            + ("" if is_venv else "\n提示：未使用项目 venv；V1 直接调用 LS-DYNA，仍可运行。")
        ),
    )


# ---------------------------------------------------------------------------
def check_ansys_dyna_core() -> None:
    section("[2] 检查 ansys.dyna.core 是否可导入")
    try:
        import ansys.dyna.core as dyna_core
    except ImportError as exc:
        warn(
            "import ansys.dyna.core",
            f"ImportError: {exc}\nPyDYNA 是可选项；当前 V1 脚本直接调用 lsdyna，不依赖它。",
        )
        return
    record(
        "import ansys.dyna.core",
        True,
        detail=f"模块路径 = {getattr(dyna_core, '__file__', '<namespace pkg>')}",
    )


# ---------------------------------------------------------------------------
def check_pydyna_version() -> None:
    section("[3] PyDYNA 版本")
    try:
        ver = pkg_version("ansys-dyna-core")
        record("ansys-dyna-core (importlib.metadata)", True, f"版本 = {ver}")
    except PackageNotFoundError:
        warn(
            "ansys-dyna-core (importlib.metadata)",
            "未找到包 ansys-dyna-core；只有需要 PyDYNA API 时才必须安装。",
        )

    try:
        import ansys.dyna.core as dyna_core

        mod_ver = getattr(dyna_core, "__version__", None)
        if mod_ver:
            record(
                "ansys.dyna.core.__version__",
                True,
                f"版本 = {mod_ver}",
            )
        else:
            record(
                "ansys.dyna.core.__version__",
                True,
                "模块未暴露 __version__（不影响使用）",
            )
    except ImportError:
        pass


# ---------------------------------------------------------------------------
def check_cursor_rule() -> None:
    section("[4] 检查 Cursor 规则文件 .cursor/rules/pydyna.mdc")
    rule_path = Path(".cursor/rules/pydyna.mdc")
    abs_path = rule_path.resolve()
    if rule_path.is_file():
        size = rule_path.stat().st_size
        record(
            ".cursor/rules/pydyna.mdc",
            True,
            f"绝对路径 = {abs_path}\n大小     = {size} bytes",
        )
    else:
        record(
            ".cursor/rules/pydyna.mdc",
            False,
            f"未找到（相对 cwd = {Path.cwd()}）\n绝对路径推断 = {abs_path}",
        )


# ---------------------------------------------------------------------------
def check_lsdyna_solver() -> None:
    section("[5] 检查 LS-DYNA 求解器（shutil.which('lsdyna')）")
    solver = shutil.which("lsdyna")
    if solver:
        record("shutil.which('lsdyna')", True, f"求解器路径 = {solver}")
    else:
        record(
            "shutil.which('lsdyna')",
            False,
            "未在 PATH 中找到 lsdyna。请检查 LS-DYNA 安装目录是否已加入 PATH。",
        )


# ---------------------------------------------------------------------------
def check_input_files() -> None:
    section("[6] 检查最小输入与运行脚本（参考项）")
    for name in ("test.k", "run_test.ps1"):
        p = Path(name)
        if p.is_file():
            record(name, True, f"大小 = {p.stat().st_size} bytes")
        else:
            record(name, False, f"未找到 -> {p.resolve()}")


# ---------------------------------------------------------------------------
def main() -> int:
    print(SEP)
    print("  LS-DYNA + PyDYNA  环境最小验证")
    print(f"  cwd = {Path.cwd()}")
    print(SEP)

    check_python()
    check_ansys_dyna_core()
    check_pydyna_version()
    check_cursor_rule()
    check_lsdyna_solver()
    check_input_files()

    section("汇总")
    n_ok = sum(1 for _, ok, _ in RESULTS if ok)
    n_total = len(RESULTS)
    print(f"  通过 {n_ok} / {n_total}")
    failed = [name for name, ok, _ in RESULTS if not ok]
    if failed:
        print("  未通过项：")
        for name in failed:
            print(f"    - {name}")
        return 1
    if WARNINGS:
        print("  提示项（不阻止 V1 运行）：")
        for name, _ in WARNINGS:
            print(f"    - {name}")
    print("  全部通过，环境就绪。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
