from swanlab.data.system import __get_command, __get_git_branch_and_commit, __get_pip_requirement

if __name__ == "__main__":
    print(__get_command())
    print(__get_git_branch_and_commit())
    print(__get_pip_requirement())
