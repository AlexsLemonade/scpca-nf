# CLAUDE.md

You are working in a Childhood Cancer Data Lab repository.
This project is Nextflow pipeline that processes data for the [Single-cell Pediatric Cancer Atlas (ScPCA) Portal](https://scpca.alexslemonade.org).

## Rules

- Do NOT run `git commit`, `git push`, or any other git write commands. The human will handle all git operations.
- Do NOT create, read, or modify `.env` files or any file containing credentials, tokens, or secrets.
- Do NOT make network requests to external services (no `curl`, `wget`, API calls, etc.). However, you may use `WebSearch` and `WebFetch` tools to retrieve documentation or code examples relevant to the task at hand.
- Do NOT install system-level packages or modify system configuration.
- Do NOT access or reference any data that is not already present in this repository.

## Expectations

- When the user asks a question about existing code or design ("why is X here?", "do we need Y?", "should this be Z?"), treat it as a discussion prompt and not as a request to make changes. Only modify code after the user has explicitly agreed on a course of action following discussion.
- When you have made changes, stop and describe what you changed and why so the human can review.
- If you are unsure whether an action is allowed, ask rather than proceeding.
- Keep changes small and incremental. Prefer multiple small steps over one large batch of changes. This helps to ensure the human can make tightly-scoped commits.

## Formatting guidelines


### Markdown and text files

- Use `#` for top-level headings, `##` for second-level, and so on. Do not skip heading levels.
- Use one sentence per line for regular text to improve readability in diffs and version control.
- Use fenced code blocks with language identifiers for code snippets (e.g., `python`, `bash`).
- Use bullet points or numbered lists for clarity when listing items or steps.
