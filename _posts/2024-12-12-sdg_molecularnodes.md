---
layout: post
title: "NumFocus Small Development Grant: Advancing Molecular Visualization with MolecularNodes"
---

<img
src="{{site.images}}/numfocus.png"
title="NumFOCUS Foundation" alt="NumFOCUS Foundation"
style="float: right; width: 10em;" />

We are thrilled to announce that MDAnalysis has been awarded a [Small Development Grant][SDG] by NumFocus to enhance scientific molecular rendering with [MolecularNodes][MN] in 2025. This initiative is a collaborative effort between [Yuxuan Zhuang][Yuxuan] and [Brady Johnston][Brady].

MolecularNodes facilitates the seamless import and visualization of structural biology data within [Blender][blender], leveraging Blender’s industry-leading visualization and animation tools. Molecular Nodes has garnered widespread excitement among scientists for its ability to create stunning and informative molecular visualizations. The current version of Molecular Nodes has an under-developed scripting interface, inhibiting the potential for automated molecular rendering. Our project aims to address these limtiations by developing a robust API, enabling users to render molecular structures with straightforward and customizable code.

<a href="https://github.com/yuxuanzhuang/ggmolvis">
    <img src="{{site.images}}/mn_example.png" 
         title="MolecularNodes" 
         alt="MolecularNodes" 
         style="display: block; margin: auto; height: 20em;" />
</a>

### Project Overview

The development will proceed in three key stages:

1. **API Development**: We will create a stable API for Molecular Nodes, empowering users to automate molecular rendering with minimal effort.
   
2. **Interactive Jupyter Integration**: A Jupyter widget will be built to integrate with MDAnalysis, providing an interactive environment for controlling and rendering molecular objects directly within notebooks via Blender.

3. **Advanced Visualization Tools**: We will develop tools for visualizing basic geometric features and even complex analysis results from MDAnalysis.

### Be Part of the Process!

We invite you to join our [Discord channel][discord] to share your ideas and feedback as we build these tools. If you'd like to be a beta user, let us know—-your input will help shape the future of Molecular Nodes! Stay tuned for updates and sneak peeks of our progress.

Thank you to [NumFocus][NumFocus] for supporting this exciting project!

- @yuxuanzhuang @bradyajohnston

[MN]: https://github.com/BradyAJohnston/MolecularNodes
[SDG]: https://numfocus.org/programs/small-development-grants
[NumFocus]: https://numfocus.org/
[Yuxuan]: https://github.com/yuxuanzhuang/
[Brady]: https://github.com/BradyAJohnston/
[discord]: https://discord.com/channels/807348386012987462/1256156074008903741
[blender]: https://www.blender.org/
