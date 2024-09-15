import os
import sys
from html import escape


def format_size(size):
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024.0:
            return f"{size:.2f} {unit}"
        size /= 1024.0


def generate_tree_html(directory):
    css = """
    .tree {
      --spacing: 1.5rem;
      --radius: 6px;
    }

    .tree li {
      display: block;
      position: relative;
      padding-left: calc(2 * var(--spacing) - var(--radius) - 2px);
    }

    .tree ul {
      margin-left: calc(var(--radius) - var(--spacing));
      padding-left: 0;
      margin-top: 0;
    }

    .tree ul li {
      border-left: 2px solid #ddd;
    }

    .tree ul li:last-child {
      border-color: transparent;
    }

    /* Horizontal lines */
    .tree ul li::before {
      content: '';
      display: block;
      position: absolute;
      top: calc(var(--spacing) / -2);
      left: -2px;
      width: calc(var(--spacing) + 2px);
      height: calc(var(--spacing) - 3px);
      border: solid #ddd;
      border-width: 0 0 2px 2px;
    }

    .tree summary {
      display: block;
      cursor: pointer;
    }

    .tree summary::marker,
    .tree summary::-webkit-details-marker {
      display: none;
    }

    .tree summary:focus {
      outline: none;
    }

    .tree summary:focus-visible {
      outline: 1px dotted #000;
    }

    /* Dots */
    .tree li::after,
    .tree summary::before {
      content: '';
      display: block;
      position: absolute;
      top: calc(var(--spacing) / 2 - var(--radius) - 2px);
      left: calc(var(--spacing) - var(--radius) - 1px);
      width: calc(2 * var(--radius));
      height: calc(2 * var(--radius));
      border-radius: 50%;
      background: #ddd;
    }

    .tree summary::before {
      z-index: 1;
      background: #696 url('expand-collapse.svg') 0 0;
    }

    .tree details[open] > summary::before {
      background-position: calc(-2 * var(--radius)) 0;
    }

    body {
      font-family: Arial, sans-serif;
      margin: 20px;
    }
    a {
      text-decoration: none;
      color: #0366d6;
    }
    a:hover {
      text-decoration: underline;
    }
    .file-size {
      color: #6a737d;
      font-size: 0.85em;
      margin-left: 0.5em;
    }
    """

    html = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="UTF-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1.0">',
        "<title>Directory Tree</title>",
        f"<style>{css}</style>",
        "</head>",
        "<body>",
        "<h1>Directory Tree</h1>",
        '<div class="tree">',
        "<li>",
        " <details open>",
        " <summary>/</summary>",
    ]

    def traverse_directory(path, depth=0):
        items = os.listdir(path)
        items.sort()

        if depth == 0:
            html.append("<ul>")

        for index, item in enumerate(items):
            full_path = os.path.join(path, item)
            relative_path = os.path.relpath(full_path, directory)

            if os.path.isdir(full_path):
                html.append("<li>")
                html.append(f"<details open><summary>{escape(item)}</summary>")
                html.append("<ul>")
                traverse_directory(full_path, depth + 1)
                html.append("</ul>")
                html.append("</details>")
                html.append("</li>")
            else:
                file_size = os.path.getsize(full_path)
                formatted_size = format_size(file_size)
                html.append(
                    f'<li><a href="{escape(relative_path)}">{escape(item)}</a>'
                    f'<span class="file-size">({formatted_size})</span></li>'
                )

        if depth == 0:
            html.append("</ul>")

    traverse_directory(directory)
    html.extend(["</li>", "</div>", "</body>", "</html>"])
    return "\n".join(html)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <directory_path>")
        sys.exit(1)

    directory = sys.argv[1]
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory")
        sys.exit(1)

    html_content = generate_tree_html(directory)
    output_file = "index.html"

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"HTML file '{output_file}' has been generated successfully.")
