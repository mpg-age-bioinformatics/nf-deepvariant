# nf-deepvariant

Running the workflow

```
RELEASE=1.0.0
PROFILE=local

nextflow run mpg-age-bioinformatics/nf-deepvariant -r ${RELEASE} -params-file params.json -entry images -profile ${PROFILE} && \
nextflow run mpg-age-bioinformatics/nf-deepvariant -r ${RELEASE} -params-file params.json -profile ${PROFILE} 
```

## Contributing

Make a commit, check the last tag, add a new one, push it and make a release:

```
git add -A . && git commit -m "<message>" && git push
git describe --abbrev=0 --tags
git tag -e -a <tag> HEAD
git push origin --tags
gh release create <tag>
```
