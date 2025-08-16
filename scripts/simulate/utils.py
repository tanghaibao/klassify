import shutil
import subprocess as sp

from pathlib import Path

from jcvi.apps.base import logger


def check_tool(name: str):
    """
    Ensure a tool is in the PATH.
    """
    if shutil.which(name) is None:
        raise RuntimeError(f"Required tool '{name}' not found in PATH")


def run(cmd, *, cwd=None, stdout=None, stdin=None):
    """
    Run a command.
    """
    logger.info("RUN: %s", " ".join(map(str, cmd)))
    sp.run(cmd, cwd=cwd, stdin=stdin, stdout=stdout, check=True)


def run_pipe(cmd1, cmd2, *, cwd=None, stdout_path: Path | None = None):
    """
    Run a pipeline of commands.
    """
    logger.info("PIPE: %s | %s", " ".join(map(str, cmd1)), " ".join(map(str, cmd2)))
    p1 = sp.Popen(cmd1, cwd=cwd, stdout=sp.PIPE)
    with open(stdout_path, "wb") as fh:
        p2 = sp.Popen(cmd2, cwd=cwd, stdin=p1.stdout, stdout=fh)
        p1.stdout.close()
        rc1 = p1.wait()
        rc2 = p2.wait()
        if rc1 != 0 or rc2 != 0:
            raise RuntimeError(f"Pipeline failed (rc1={rc1}, rc2={rc2})")
