# RepeatViewer User Manual

Read it as web pages (`index.html` in this folder, also under Help in the app) or as a [PDF](RepeatViewer_User_Manual.pdf). The chapters below are the sources.

1. [Introduction](01-introduction.md)
2. [Installation](02-installation.md)
3. [Quick start: one curation pass](03-quick-start.md)
4. [The interface](04-interface.md)
5. [Sequences](05-sequences.md)
6. [Annotations](06-annotations.md)
7. [Search](07-search.md)
8. [Alignments](08-alignments.md)
9. [Mapping](09-mapping.md)
10. [Trees](10-trees.md)
11. [Repeat analysis](11-repeat-analysis.md)
12. [Primers and PCR](12-primers.md)
13. [Automation](13-automation.md)
14. [Settings](14-settings.md)
15. [Reference](15-reference.md)

To change the manual, edit the chapter files, then rebuild:

```
python docs/manual/tools/build_manual.py --pdf
```

After a change to a window, take its screenshot again:

```
python docs/manual/tools/take_screenshots.py --only NAME
```
