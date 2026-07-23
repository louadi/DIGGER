# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

DIGGER 2.0 is a Django web application for exploring how alternative splicing affects protein–protein interactions (PPIs). It joins PPIs and domain–domain interactions (DDIs) into a single graph and maps interacting residues down to exons. Public instance: https://exbio.wzw.tum.de/digger

The repo has two independent halves:
- `container/` — the Django web app (served in production via Docker).
- `preprocess/` — an **offline** pipeline that predicts DDIs for a new organism and generates the data files the web app consumes. It is only run manually when adding an organism; it is not part of the running service.

## Commands

All web-app commands run from `container/` (or via `docker-compose exec web ...`).

- **Run the app locally (SQLite, no Docker):** `DJANGO_SETTINGS_MODULE=DomainExplorer.settings.development python manage.py runserver`. `manage.py` defaults to the `development` settings module.
- **Tests:** `pytest` from `container/`. Add `-s` for prints, `-n auto` for parallel (pytest-xdist). Tests need `DJANGO_SETTINGS_MODULE` set (there is no conftest/ini entry for it) — use `development` for SQLite. `domain/tests.py` is a heavily parametrized sweep of NEASE runs over the fixtures in `domain/test_files/` and uses `@pytest.mark.django_db`.
- **Load the biological data into the DB:** `python manage.py import_datasets` (bulk-loads CSVs — see data model below).
- **Standard Django:** `python manage.py migrate`, `collectstatic --no-input`, `createsuperuser`.

Full first-time deploy sequence (download data from Zenodo, `docker-compose up -d --build`, migrate, import, collectstatic, createsuperuser) is in `README.md`. Data archives (`data.zip`, `nease_data.zip`) unzip into `container/domain/` and `container/domain/nease/` respectively — they are not in git.

Extending with a new organism's DDIs: `conda env create -f DIGGER_env.yml`, then `cd preprocess && python main_ddi_extend.py`, driven by `preprocess/sourcedata/database_sources.yml` (Linux-only; requires several source files listed in `main_ddi_extend.py`).

## Architecture

### Django project layout
- `DomainExplorer/` is the Django project; `domain/` is the single app holding all logic.
- **Settings are split and layered:** `common.py` + `i18n.py` are shared; `production.py` (Postgres, `ALLOWED_HOSTS=['*']`) and `development.py` (SQLite at `run/dev.sqlite3`) each import them. `settings/__init__.py` re-exports `production`, so `import DomainExplorer.settings` = production. Docker sets `DJANGO_SETTINGS_MODULE=DomainExplorer.settings.production`; `manage.py` defaults to `development`.
- Every settings module defines `DATABASE_ENGINE`, a **SQLAlchemy engine** used *alongside* the Django ORM (see data model).
- Stack is intentionally old: Python 3.7, Django 2.2. Match existing style.

### Two parallel data stores (this is the key thing to understand)
Biological data does **not** live in the Django ORM. There are only three ORM models (`Gene`, `Domain`, `NeaseSavedRun` in `domain/models.py`). The real data is served two ways, and **both must be populated per organism**:

1. **In-memory pickle graphs**, loaded once at import time by `domain/Process/load_data.py`. It scans `domain/data/`, and for each organism folder builds entries in module-level dicts (`PPI_all`, `g2d_all`, `protein_all`, `DDI_all`, `DomainG_all`, `gid2name_all`). Used for graph/network operations (networkx).
2. **SQL tables**, populated by `manage.py import_datasets`, which bulk-loads per-organism CSVs (`PPI_interface_mapped_to_exon.csv`, `final.csv`, `gene_info.csv`, `gene_name2entrez_id.csv`, `Pfam-A.clans.tsv`) into tables named `ppi_data_<trivial>`, `exons_to_domains_data_<trivial>`, etc. The `Process/` modules query these with `pd.read_sql_query(..., con=settings.DATABASE_ENGINE)` — raw SQL through SQLAlchemy, not the ORM.

### Organism naming convention
Organism data folders are named `Scientific name[trivial]`, e.g. `Homo sapiens[human]`, `Mus musculus[mouse]`. The bracketed **trivial name is the key used everywhere** (dict keys, SQL table suffixes, URL segments). It is extracted with `organism.split("[")[1][:-1]`. Adding an organism = drop a correctly-named folder with all expected files into `domain/data/`, then re-run `import_datasets` and restart (to re-trigger the import-time pickle load).

### Request flow
- `DomainExplorer/urls.py` mounts static pages (`DomainExplorer/views.py`, class-based) and includes `domain/urls.py`.
- `domain/urls.py` routes organism-scoped URLs (organism is a path segment) to function-based views in `domain/views.py`: `gene`, `exon`, `transcript`, `isoform_level`, `exon_level`, `network`, `Multi_proteins`, `setup_nease`, `multiple_queries`. Legacy DIGGER 1.0 URLs 301-redirect to their `human` equivalents.
- Views are thin; computation lives in `domain/Process/` (`exon.py`, `gene.py`, `transcript.py`, `exonstodomain.py`, `network_analysis.py`, `mutliple_query.py` [note the misspelling — keep it], `proteininfo.py`, `process_data.py`, `nease_output.py`).
- Admin is served at `secure-admin/`, not `admin/`.

### NEASE
`domain/nease/` is a vendored copy of the NEASE library (enrichment analysis of splicing events), entered via `from .nease import run` (`domain/Process/nease_output.py` wraps it). NEASE run results are persisted as files and tracked by the `NeaseSavedRun` model, which records a retention bucket (`saved_for_days`) and computes `days_left()`.

### Deployment specifics
- `docker-compose.yml`: `web` (gunicorn via `docker/web/start.sh`), `db` (postgres:12, exposed on host port **5461**), `nginx` (published on host port **8061**).
- **Subpath hosting:** `NGINX_PUBLISHED_PATH` (from `.env`) prefixes `STATIC_URL`, `MEDIA_URL`, and the session/CSRF cookie paths, so multiple Django instances can share a host. `BASE_URL` must be set correctly or reference links break.
- A daily cron (`docker/web/auto_clean.sh`) deletes aged NEASE result files from `nease_events/{zero,seven,thirtyone,onehundredeightysix}_days` folders by mtime, plus old files in `run/media`.
- `run/SECRET.key` is auto-generated on first start if absent.
