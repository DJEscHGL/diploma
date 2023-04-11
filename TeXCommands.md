# Кастомные команды для TeX

### Главы без нумерации

```tex
\anonsection

\newcommand{\anonsection}[1]{\section*{#1}\addcontentsline{toc}{section}{#1}}
```

![](Files/Assets/Screenshots/CC1.png)
