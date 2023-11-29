"""Utilities for getting information from git SCM."""
import git


def get_git_commit_hash() -> str:
    """Get the git commit hash of the repository in which this code is found.

    Returns
    -------
    str:
        Commit hash of the current repo. Returns "unknown" if the code is not
        in a git repository. Appends "(dirty)" if there are uncommitted
        changes.

    """
    try:
        repo = git.Repo(__file__, search_parent_directories=True)
    except git.InvalidGitRepositoryError:
        return "unknown"

    commit_hash = str(repo.head.commit)
    if repo.is_dirty():
        commit_hash += " (dirty)"

    return commit_hash


def simplify_remote(remote: str) -> str:
    """Simplify a remote URL if necessary to keep it below 64 characters."""
    if len(remote) <= 64:
        return remote

    # Strip from start and end
    if remote.endswith(".git"):
        remote = remote[:-4]
    for start_str in ["http://", "https://", "git@"]:
        if remote.startswith(start_str):
            remote = remote[len(start_str):]

    if len(remote) <= 64:
        return remote

    remote = remote.replace("github.com:", "github:")
    remote = remote.replace("github.com/", "github:")

    if len(remote) > 64:
        raise ValueError(
            "Cannot simplify URL of the remote to be 64 characters or fewer."
        )

    return remote


def get_git_remote_url(simplify: bool = True) -> str:
    """Get the remote URL of the git repository in which this code is found.

    Returns
    -------
    str:
        Remote URL of the git repository. Returns "unknown" if the code is not
        in a git repository or the repository has no remotes. If there are
        multiple remotes, the first is used.

    """
    try:
        repo = git.Repo(__file__, search_parent_directories=True)
    except git.InvalidGitRepositoryError:
        return "unknown"

    if len(repo.remotes) == 0:
        return "unknown"

    remote = repo.remotes[0].url

    # Convert an ssh remote to a https remote
    if remote.startswith("git@"):
        remote = (
            remote
            .replace(":", "/")
            .replace("git@", "https://")
        )

    if simplify:
        remote = simplify_remote(remote)

    return remote
