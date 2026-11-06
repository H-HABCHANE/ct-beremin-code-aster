"""Lance Code_Aster dans son image Docker, sans dependance sur l'hote."""
from __future__ import annotations

import os
import shutil
import subprocess
import uuid
from pathlib import Path

IMAGE = "simvia/code_aster:17.4.0"


def write_export(case_dir: Path, comm: str, data: dict[int, str], results: dict[int, str], memory_mb: int = 3000, ncpus: int = 1) -> Path:
    lines = [
        "P time_limit 1800",
        f"P memory_limit {memory_mb}",
        f"P ncpus {ncpus}",
        "P mpi_nbcpu 1",
        f"F comm {comm} D 1",
    ]
    lines += [f"F libr {name} D {unit}" for unit, name in data.items()]
    lines += [f"F libr {name} R {unit}" for unit, name in results.items()]
    lines.append("F mess message.mess R 6")
    export = case_dir / "calcul.export"
    export.write_text("\n".join(lines) + "\n")
    return export


def run(case_dir: Path, timeout_s: float | None = None) -> None:
    """Lance le calcul. En cas de depassement, le conteneur est arrete : il ne survit pas au script."""
    name = f"aster-{uuid.uuid4().hex[:12]}"
    cmd = [
        "docker", "run", "--rm", "--name", name,
        "--user", f"{os.getuid()}:{os.getgid()}",
        "--env", "HOME=/tmp",
        "--mount", f"type=bind,source={case_dir.resolve()},target=/calc",
        "--workdir", "/calc",
        IMAGE, "bash", "-c", "source /opt/activate.sh && run_aster calcul.export",
    ]
    log = case_dir / "docker.log"
    with log.open("w") as handle:
        try:
            proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, text=True, timeout=timeout_s)
        except subprocess.TimeoutExpired:
            subprocess.run(["docker", "kill", name], capture_output=True)
            raise RuntimeError(f"Code_Aster : duree limite de {timeout_s} s depassee, voir {log}") from None
    if proc.returncode != 0:
        raise RuntimeError(f"Code_Aster a echoue (code {proc.returncode}), voir {case_dir}/message.mess")


def copy_comm(source: Path, case_dir: Path) -> str:
    shutil.copy(source, case_dir / source.name)
    return source.name
