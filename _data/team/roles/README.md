# Roles

This directory contains YAML files, one per role.

## Editing

Each YAML file should contain the details of one role only.
Each role is automatically picked up by the template. To add or remove a role, simply add or remove a file. Do not add non-role documents to this directory that can be parsed by Jekyll.

Each file *must* contain the following tags at top level:

* `role`: The name of the role (e.g. Core library maintenance)

Each role can optionally have the following tags:
* `description`: text describing the role. This can be written in Markdown.
* `tasks`: a list of strings describing tasks. This is rendered as a list.

Each role can also optionally have lead/member tags, if there are no subroles. If there are subroles, the below are ignored:
* `current_leads`: a list of names who are current leads
* `current_members`: a list of names who are current members in the subgroup
* `historical_leads`: a list of names who were previously leads
* `historical_members`: a list of names who were previously members in the subgroup

Onto subroles: each role can optionally define `subroles`. The `subroles` tag should be a **list** where each item *must have*:

* `subrole`: The name of the subrole (e.g. Issue management)

Each subrole can *optionally have*:
* `current_leads`: a list of names who are current leads
* `current_members`: a list of names who are current members in the subgroup
* `historical_leads`: a list of names who were previously leads
* `historical_members`: a list of names who were previously members in the subgroup
* `description`: text describing the role. This can be written in Markdown.
* `tasks`: a list of strings describing tasks. This is rendered as a list.


## Parsing

The files here are read first by ``_includes/team_table.html`` and ``_includes/roles_description.html``.
Jekyll automatically reads YAML, JSON, and other formatted files in the ``_data`` directory
as objects that can be interacted with using the Liquid templating language.
For example, the ``name`` variable of the ``code_of_conduct.yml`` file is accessible as ``site.data.team.roles.code_of_conduct.name`` in HTML.

