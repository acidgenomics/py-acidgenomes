from pathlib import Path

import pytest

import acidgenomes._genome_version as genome_version
import acidgenomes._mapping as mapping
import acidgenomes.download._ensembl as ensembl_download


def test_current_ensembl_version_uses_version_endpoint(monkeypatch):
    calls: list[str] = []

    def mock_fetch_text(url: str) -> str:
        calls.append(url)
        return "116\n"

    monkeypatch.setattr(genome_version, "fetch_text", mock_fetch_text)

    assert genome_version.current_ensembl_version() == 116
    assert calls == ["https://ftp.ensembl.org/pub/VERSION"]


def test_current_ensembl_version_rejects_non_version_response(monkeypatch):
    monkeypatch.setattr(genome_version, "fetch_text", lambda url: "Ensembl Release 116")

    with pytest.raises(RuntimeError, match="Failed to parse Ensembl release from VERSION"):
        genome_version.current_ensembl_version()


def test_download_ensembl_genome_uses_https_release_urls(monkeypatch, tmp_path: Path):
    urls: list[str] = []

    def mock_download(url: str, dest: Path) -> Path:
        urls.append(url)
        return dest

    monkeypatch.setattr(ensembl_download, "_download_file", mock_download)

    ensembl_download.download_ensembl_genome(
        organism="Homo sapiens",
        genome_build="GRCh38",
        release=116,
        output_dir=tmp_path,
    )

    assert len(urls) == 4
    assert all(url.startswith("https://ftp.ensembl.org/pub/release-116/") for url in urls)


def test_map_ensembl_release_to_url_uses_rest_for_current():
    assert mapping.map_ensembl_release_to_url() == "https://rest.ensembl.org"


def test_map_ensembl_release_to_url_uses_canonical_archive_index(monkeypatch):
    calls: list[str] = []

    def mock_fetch_text(url: str) -> str:
        calls.append(url)
        return "Ensembl 100 Apr 2020"

    monkeypatch.setattr(mapping, "fetch_text", mock_fetch_text)

    assert mapping.map_ensembl_release_to_url(100) == "https://apr2020.archive.ensembl.org"
    assert calls == ["https://www.ensembl.org/info/website/archives/index.html"]
