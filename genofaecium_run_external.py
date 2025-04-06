import subprocess


def run_shell_command_tool_check(command):
    """Run a shell command and return success flag and output."""
    try:
        result = subprocess.run(command, shell=True, executable="/bin/bash",
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        output = result.stdout.strip() or result.stderr.strip()
        success = result.returncode == 0
        return success, output
    except Exception as e:
        return False, str(e)


def run_tool_check_in_env(env_path, tool_command):
    """Run a tool in a conda-packed environment."""
    bash_cmd = f"""
    source {env_path}/bin/activate
    {tool_command}
    source ~/.bashrc
    """
    return run_shell_command_tool_check(bash_cmd)



def run_shell_command_tool_job(command):
    """Run a shell command and return success flag and output."""
    try:
        result = subprocess.run(command, shell=True, executable="/bin/bash",
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        output = ''
        success = result.returncode == 0
        return success, output
    except Exception as e:
        return False, str(e)


def run_tool_job_in_env(env_path, tool_command):
    """Run a tool in a conda-packed environment."""
    bash_cmd = f"""
    source {env_path}/bin/activate
    {tool_command}
    source ~/.bashrc
    """
    print("[GenoFaecium] [genofaecium_run_external] executing the following command:  " + bash_cmd)
    return run_shell_command_tool_job(bash_cmd)


