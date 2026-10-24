# Task front end: "Choose Your Friends Wisely!"

The participant-facing code for this study's task, as it runs live on
[MySocialBrain.org](https://mysocialbrain.org), an ongoing citizen-science
platform hosting many unrelated studies. This is shared **for transparency
about exactly how the task was presented and how responses were
recorded** -- it is not a standalone, runnable copy of the task, and no
attempt has been made to make it one.

## What's here

- **`net_info.html`** -- the study's public landing page (what a visitor
  sees before starting: study description, "Begin" button).
- **`pages/net_consent.html`** -- the IRB consent form.
- **`pages/net_instruct.html`** -- task instructions, shown once per
  participant before the first trial.
- **`pages/net_test.html`** -- the trial page itself: the network
  visualization area and the "Confirm" button, repeated once per trial.
- **`pages/net_results.html`** -- the post-task debrief/results page
  (explains the study's purpose and hypotheses, shows the participant's
  own results).
- **`js/network_task.js`** -- the task's client-side logic, extracted
  verbatim from the platform's shared `js/studylogic.js` (which holds the
  equivalent code for every study on the platform in one file). Renders
  each trial's stimulus graph with D3, records the participant's choice
  and response time, and (after the final trial) renders the per-subject
  feedback visualization.
- **`css/network_task.css`** -- the node/link styling this task's
  visualization uses, extracted from the platform's shared `css/styles.css`.
- **`img/network.png`** -- the example-network screenshot shown on the
  instructions and landing pages.

## What's deliberately NOT here, and why

This is genuinely not runnable outside the live platform, by design:

- **No server-side code.** Everything that generates stimulus graphs,
  serves them to the client, or records responses (`get_networks.php`,
  `netresponses.php`, `getvis.php`, `getindivis.php`, and the database
  layer behind them) lives server-side and is not part of this archive.
  The stimulus graphs themselves -- their generation, not their serving --
  are in `generate_networks.R` in the parent repository.
- **No shared platform framework.** The task code above calls into
  `window.EM` (`modules.post`, `modules.getPage`, `modules.randInt`,
  session/auth handling, page routing, etc.) -- the platform's shared
  JavaScript infrastructure used by every study on the site, not specific
  to this one. That framework isn't included, so `network_task.js` won't
  run as-is.
- **No shared page chrome.** The site-wide header, footer, navigation,
  and the generic libraries the platform loads on every page (jQuery,
  Bootstrap, D3 itself, etc.) aren't included either.

## Provenance

Copied from the live site's front-end source (`D:\Dropbox\webprojects\
mysocialbrain` at the time of this archive), touching only HTML, CSS,
JS, and image files -- no PHP or SQL/database-layer files were read or
copied in preparing this folder.
