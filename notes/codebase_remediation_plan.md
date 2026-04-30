# Codebase Remediation Plan

A working document combining lessons from a talk on AI-era software fundamentals with a detailed remediation strategy for a months-old AI-developed codebase running on an airgapped HPC.

---

## Part 1: Notes from the Talk — "Software Fundamentals Matter More Than Ever"

### Core Thesis

Good codebases matter more than ever in the AI era. Bad code is now the most expensive it's ever been because AI thrives in well-structured codebases but produces garbage in poorly-structured ones. The "specs-to-code" movement (write a spec, generate code, ignore the code, iterate the spec) is essentially vibe coding rebranded — and it produces progressively worse code with each iteration.

### The Six Failure Modes & Solutions

**1. The AI didn't do what I wanted**

The problem is that you and the AI don't share a "design concept" (Frederick P. Brooks' term for the invisible, ephemeral idea of what you're building). The fix is a **"Grill Me" skill**: instruct the AI to interview you relentlessly about every aspect of the plan, walking down each branch of the design tree, until you reach shared understanding. This can mean 40–100+ questions before the AI is satisfied. The speaker prefers this over Claude Code's default plan mode, which is too eager to create an asset rather than reach genuine alignment first.

**2. The AI is way too verbose**

This is a language gap, similar to working with a domain expert in an unfamiliar field. The solution comes from Domain-Driven Design: establish a *ubiquitous language* — a shared vocabulary used in code, conversations, and AI prompts. The speaker built a skill that scans the codebase, extracts terminology, and generates a markdown file of terms. Keeping this open during planning improves AI thinking traces, reduces verbosity, and aligns implementation with intent.

**3. The AI built the right thing, but it doesn't work**

Feedback loops are essential — TypeScript, browser access for frontend work, automated tests. But LLMs misuse feedback loops by "outrunning their headlights" (Pragmatic Programmer): producing huge amounts of code before checking anything. The rate of feedback is your speed limit. TDD forces small, deliberate steps: write a test, make it pass, refactor.

**4. Testing is hard**

Testing requires interdependent decisions (unit size, what to mock, which behaviors to test). Good codebases are testable codebases. The fix is John Ousterhout's *deep modules*: lots of functionality hidden behind simple interfaces, rather than many shallow modules with complex interfaces. Shallow-module codebases are exactly what AI tends to produce, and they're hard for AI to navigate. The speaker has an "improve codebase architecture" skill that finds related code and wraps it in deep modules — creating simple boundaries that reward TDD.

**5. Your brain can't keep up**

Even with working feedback loops, shipping more code than ever is exhausting. Deep modules let you treat implementations as *gray boxes*: you design the interface carefully but delegate the internals to the AI (with caveats for critical areas like finance). You test from the outside and verify behavior, saving cognitive load.

**6. Staying aware of structure**

Every plan, PRD, and conversation should reference specific modules and interfaces being modified. Per Kent Beck: "invest in the design of the system every day." Specs-to-code does the opposite — it divests from design.

### The Big Takeaway

Think of AI as a brilliant tactical programmer — a sergeant making changes on the ground. You're the strategist above it. That strategic role requires the same software fundamentals developers have used for 20+ years: design thinking, deep modules, ubiquitous language, TDD, and continuous investment in architecture.

### Referenced Resources

- *A Philosophy of Software Design* by John Ousterhout
- *The Pragmatic Programmer*
- *The Design of Design* by Frederick P. Brooks
- Domain-Driven Design (ubiquitous language concept)
- Kent Beck on continuous design investment
- Speaker's skills repo (GitHub) — search for "mac PCO skills" or similar
- aihero.dev (speaker's newsletter)

---

## Part 2: Diagnosed Problems in the Codebase

### Problems Already Identified

1. **Over-complex code** resulting from continual fixes, rather than writing simple correct code to achieve the core task.
2. **Repetition of functionality** across many scripts — several functions may do the same thing in different ways.
3. **Dead code** remaining from previous iterations that is never used in the existing codebase.
4. **Over-verbose comments** describing every change over many iterations.
5. **Inconsistent file sizes** — some scripts are 6,000 lines, others 50 lines, with no coherent grouping principle.
6. **Terminology and naming inconsistencies** across files; many variables with similar names that are unclear to human developers.
7. **Complete lack of unit tests.**

### Other Problems Likely Present

Beyond the seven above, long-running AI-developed codebases typically also exhibit:

**Architectural drift.** Functions and modules that no longer reflect the original mental model. The codebase has accumulated layers like sedimentary rock — early decisions buried under later ones that contradict them.

**Inconsistent error handling.** Some functions raise exceptions, some return None, some print warnings, some silently fail. AI tends to pick whatever pattern is locally convenient rather than what's globally consistent.

**Phantom abstractions.** Classes, wrappers, or helper functions created "just in case" or to solve a problem that no longer exists. They add cognitive load without earning their keep.

**Configuration sprawl.** Hardcoded values, magic numbers, and parameters scattered across files instead of centralized. Often the same constant defined in three places with slightly different values.

**Implicit coupling.** Modules that appear independent but secretly depend on each other through shared global state, file paths, environment variables, or assumptions about call order.

**Inconsistent data structures.** The same conceptual entity (e.g., a "result" or a "record") represented as a dict in one place, a dataclass in another, a tuple in a third. Conversions happening implicitly.

**Logging and debugging cruft.** Print statements, commented-out logging, debug flags that no longer work — relics of past debugging sessions.

**Untracked invariants.** Assumptions the code relies on but never checks or documents (e.g., "this list is always sorted," "this dict always has key X").

**Documentation rot.** Docstrings that describe what the function used to do, not what it does now. README files that reference deleted modules.

**Dependency bloat.** Imports of libraries used by one tiny function, or three libraries that do the same thing.

---

## Part 3: Strategic Principles

A few principles that should guide everything:

**Don't refactor and add features at the same time.** This is the cardinal rule. Functionality is the safety net. Since current functionality appears correct, preserve that absolutely. Every refactoring step should be behavior-preserving.

**A safety net is needed before refactoring can begin.** This is the chicken-and-egg problem the talk doesn't fully address. Refactoring without tests is unsafe, but writing tests for bad code is painful. The answer is *characterization tests* (also called golden-master or approval tests): tests that capture current behavior, even buggy behavior, so deviations are detectable.

**The HPC constraint changes the workflow significantly.** Without rapid Claude Code iteration on the live codebase, the approach must be more deliberate and planning-heavy — which fits the talk's philosophy well. Treat it as forced discipline.

**Work in small, reversible commits.** Even if version control isn't being used religiously now, start. Every refactoring step should be a separate commit so issues can be bisected.

---

## Part 4: HPC Workflow Setup

Given the airgapped HPC constraint, the workflow should be:

- Develop changes locally on the Mac with Claude Code (or Claude desktop/web).
- For each refactoring task: develop the change locally, copy modified files to the HPC, run tests/the actual code there, copy results back, iterate.
- Maximize what each round-trip accomplishes: do thorough planning before any code change, batch related changes, and lean heavily on static analysis (which doesn't need execution) for the early phases.

**Practical setup:**

- Get the codebase into git if it isn't already.
- Set up a way to run the code's test suite (even a trivial smoke test) on the HPC.
- Create a sync script (rsync or controlled-transfer equivalent) for the round-trips.
- Factor transfer latency into planning if the airgap requires physical media or controlled transfer.

---

## Part 5: The Remediation Plan

The early phases are diagnostic and require no code execution, which suits the HPC constraint.

### Phase 1: Map the Territory (Local, No Execution Needed)

Before changing anything, build a mental model of what exists. The **ubiquitous language skill** earns its keep first.

**Step 1.1 — Structural inventory.** Have Claude walk the codebase and produce:
- A list of all modules/scripts with line counts.
- A list of all top-level functions and classes with one-line summaries.
- A dependency graph (which modules import which).
- A rough categorization of what each module is responsible for.

This is purely a read-only operation.

**Step 1.2 — Run the ubiquitous language skill.** This will surface naming inconsistencies (problem #6) immediately. The same concept likely has three names, or one name means three things. Don't fix anything yet — just produce the glossary.

**Step 1.3 — Dead code pass with static analysis.**
- `vulture` (Python) finds unused functions and imports.
- `pyflakes` and `ruff` find unused variables.

These run locally and identify problem #3 essentially for free. Produce a list, but don't delete anything yet.

**Step 1.4 — Complexity report.** Ask Claude (and use tools like `radon` for cyclomatic complexity) to identify, for each module: the longest functions, the most deeply nested code, the most complex control flow. This maps problem #1.

**Phase 1 deliverables:** module map, glossary, dead-code list, complexity report, and — crucially — a written description of what the codebase actually does at the level of major functional areas. That last document is the *design concept*, made explicit for the first time.

### Phase 2: Establish Behavioral Tests Before Refactoring

This addresses problem #7, but more importantly creates the safety net. Don't try to write comprehensive unit tests yet — that's a Phase 5 activity. For now, write characterization tests at the highest level.

**Step 2.1 — Identify critical outputs.** Pick the 5–15 things the codebase produces that actually matter — specific computational results, output files, plots, or whatever the domain produces.

**Step 2.2 — Capture golden-master outputs.** Run the codebase on representative inputs (this happens on the HPC) and save the outputs. Treat these as the golden master.

**Step 2.3 — Write comparison scripts.** Minimal scripts that re-run those computations and compare against the golden master with appropriate tolerances (numerical comparisons need this). These are end-to-end characterization tests. They're ugly, slow, and not what you'd write in a clean codebase — but they tell you if a refactor broke anything that matters.

**Step 2.4 — Add instrumentation if needed.** If parts of the codebase don't produce easily-comparable outputs, add small instrumentation to capture intermediate state for comparison.

This phase requires HPC access since the code needs to actually run.

### Phase 3: Mechanical Cleanup (Low Risk, High Value)

Now code changes begin, but only in mechanically-safe ways. Each of these can be its own commit; characterization tests verify nothing broke.

**Step 3.1 — Delete dead code.** Use the Phase 1 list. Verify with static analysis that nothing references it. Run characterization tests.

**Step 3.2 — Strip over-verbose comments (problem #4).** Have Claude review each module and rewrite comments to describe *current* behavior only, removing change history (which belongs in git anyway). Purely a readability change.

**Step 3.3 — Fix naming inconsistencies (problem #6).** Use the glossary. Pick canonical names for each concept and use find/replace systematically. Be careful with names that are too generic (`data`, `result`) — these need context-aware renaming.

**Step 3.4 — Centralize configuration.** Pull magic numbers and paths into a config module or YAML file.

**Step 3.5 — Consolidate duplicated functionality (problem #2).** For each cluster of similar functions, pick the best implementation, replace calls to the others, then delete the rest. The ubiquitous language work makes finding these clusters much easier.

After this phase the codebase should be noticeably smaller and cleaner without any structural changes yet.

### Phase 4: Structural Refactoring Toward Deep Modules

This is where the talk's main advice applies most directly. Now problems #1 and #5 get addressed: re-organize code into deep modules with clear interfaces.

**Step 4.1 — Design the target architecture.** Use the **grill-me skill** in an unusual way: have Claude grill you on what the *ideal* architecture should look like, given what the code actually does. You're designing the target state. Produce a document describing the modules you want, their interfaces, and their responsibilities.

**Step 4.2 — Migrate one module at a time.** Use the **improve codebase architecture skill** the speaker mentions, or do it manually:
- Identify groups of related functionality currently scattered across files.
- Define the interface (the small set of functions/classes the rest of the code will call).
- Migrate code into the new module behind that interface.
- Each module migration is its own commit; characterization tests verify behavior.

**Step 4.3 — Address giant files thoughtfully.** Don't blindly split by line count. Split by responsibility, identified through the module design. Some 6,000-line files might be legitimately one module if they're cohesive; most won't be.

This phase is the most demanding round-trip-wise. Plan thoroughly with Claude locally, make changes in one batch per module, test on HPC, iterate. Resist the temptation to do many modules at once.

### Phase 5: Add Real Unit Tests

Once deep modules with clean interfaces exist, unit testing becomes tractable — exactly the talk's point. Now TDD-style discipline applies for any *future* changes, and unit tests can be backfilled for the most critical or most-changed modules.

Don't aim for 100% coverage. Aim for:
- Tests on the public interface of each deep module.
- Tests on any logic that's tricky enough to be worth verifying independently.

Characterization tests stay in place as integration-level safety nets.

### Phase 6: Establish Going-Forward Discipline

The whole point is to not end up here again. Practical commitments:

- Use the grill-me skill before any non-trivial change.
- Keep the ubiquitous language doc updated as the canonical glossary.
- Run static analysis before commits.
- Write tests for new functionality (TDD where feasible).
- Treat any module that grows past some threshold (say 500 lines) as a refactoring trigger.

---

## Part 6: Concrete First Steps for This Week

1. Get the codebase into git locally if it isn't already.
2. Set up rsync or transfer mechanism between Mac and HPC.
3. Run the structural inventory step from Phase 1 — have Claude read the codebase and produce the four documents. No execution needed; gives enormous insight before committing to any plan.
4. Run the ubiquitous language skill to produce the glossary.
5. Run static analysis tools (`ruff`, `vulture`, `radon`) for dead code and complexity reports.
6. Read all outputs carefully. They will probably reframe remediation priorities — the real problems may be different from the top-of-mind list.
7. *Then* use the grill-me skill to plan the actual refactoring, armed with concrete data about the codebase rather than impressions.
8. Only after that planning, start touching code — and the first code changes should be Phase 3 mechanical cleanups, not architectural restructuring.

---

## Part 7: Honest Scope Warning

A months-long AI-grown codebase can take weeks of disciplined remediation to clean up properly. The temptation will be to rush to architectural fixes because they feel most satisfying, but the boring early phases — mapping, characterization tests, mechanical cleanup — are what make the architectural work safe. Skipping them turns a working-but-messy codebase into a broken-and-still-messy one.

---

## Quick Reference: Skills to Use

| Skill | When to Use | Phase |
|---|---|---|
| Ubiquitous Language | Build the glossary; resolve naming inconsistencies | 1, 3 |
| Grill Me | Plan target architecture; plan any non-trivial change | 4, 6 |
| Improve Codebase Architecture | Migrate code into deep modules | 4 |

## Quick Reference: Static Analysis Tools (Python)

| Tool | Purpose |
|---|---|
| `ruff` | Fast linting; unused variables, style issues |
| `vulture` | Dead code detection (unused functions, classes, imports) |
| `pyflakes` | Lightweight error checking |
| `radon` | Cyclomatic complexity metrics |
| `mypy` | Static type checking (if type hints present) |
