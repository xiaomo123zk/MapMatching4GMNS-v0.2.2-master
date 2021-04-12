import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="MapMatching4GMNS",
    version="0.2.2",
    author="Xuesong (Simon) Zhou, Kai (Frank) Zhang, Jiawei Lu",
    author_email="xzhou74@asu.edu, zhangk2019@seu.edu.cn, jiaweil9@asu.edu",
    description="An open-source, cross-platform, lightweight, and fast Python\
                 MapMatching4GMNS engine for mapping GPS traces to the underlying network\
                 using General Modeling Network Specification (GMNS).\
                 Its most likely path finding algorithm takes about 0.02 seconds to process one GPS trace\
                 with 50 location points in a large-scale network with 10K nodes.",
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/asu-trans-ai-lab/MapMatching4GMNS",
    packages=['MapMatching4GMNS'],
    package_dir={'MapMatching4GMNS': 'MapMatching4GMNS'},
    package_data={'MapMatching4GMNS': ['bin/*']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
)
