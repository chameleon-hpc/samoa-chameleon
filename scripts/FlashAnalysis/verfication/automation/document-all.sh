#!/bin/bash

source './config.sh'

if [ "$1" == 'html' ]; then
    echo "HTML mode"
elif [ "$1" == 'pdf' ]; then
    echo "PDF mode"
else
    echo "Unknown mode"
    exit 1
fi

echo 'Starting documentation job'

rm -rf '../documentation_tmp'
for RENDER in renders/*; do
    RENDERE=$(basename -- "$RENDER")
    RENDERE=${RENDERE%.*}
    OUTPUT="../documentation_tmp/$RENDERE/$RENDERE.md"
    mkdir -p "../documentation_tmp/$RENDERE"
    echo "Writing docs for $RENDERE"
    cat "docs/$RENDERE.md" > $OUTPUT
    COMPILE="## Compilation and Execution"
    for RUN in runs/$RENDERE*; do
        RUNE=$(basename -- "$RUN")
        RUNE=${RUNE%.*}
        echo "Adding compilation from $RUNE"
        source $RUN
        source "builds/$BINARY.sh"
        for LIMITER in $LIMITERS; do
            COMPILE="$COMPILE\n\n### $DESCRIPTION using the $LIMITER limiter\n\nCompilation:\n\n\`\`\`\nswe_scenario=$BINARY\nlimiter=$LIMITER\n$BUILDFLAGS\n\`\`\`\n\nExecution:\n\n\`\`\`\n${PARAMS//\//\\/} -xdmfoutput -output_dir output\/verify\/${RUNE}_$LIMITER\n\`\`\`"
        done
    done
    sed -i "s/## Compilation and Execution/$COMPILE/g;" $OUTPUT
    echo "Adding rendering from $RENDERE"
    for LIMITER in $LIMITERS; do
        echo -e "\n### Using the $LIMITER limiter\n\n\`\`\`" >> $OUTPUT
        cat $RENDER | grep '.py' | sed "s/\\\$OUTPUT[0-9]*/results\/${RENDERE}_$LIMITER/g; s/\\\$1/$LIMITER/g;" >> $OUTPUT
        echo -e "\`\`\`" >> $OUTPUT
    done
    echo -e "\n## Results" >> $OUTPUT
    RESBASE="../results/$RENDERE/"
    for RESULT in `find "$RESBASE" -type f -name "*.svg"`; do
        NAME=${RESULT/$RESBASE/}
        echo "Adding result $NAME"
        if [ "$1" == 'pdf' ]; then
            echo -e "\n### $NAME\n\n![](${INDICES[$RENDERE]}_${RENDERE}/$NAME)\n" >> $OUTPUT
        else
            echo -e "\n### $NAME\n\n![]($NAME)\n" >> $OUTPUT
        fi
    done
done

echo 'Copying results'
cp -rf '../results/.' '../documentation_tmp'
echo 'Indexing renderings'
for RENDER in `find '../documentation_tmp' -maxdepth 1 -mindepth 1 -type d`; do
    RENDERE=$(basename -- "$RENDER")
    mv "$RENDER" "../documentation_tmp/${INDICES[$RENDERE]}_${RENDERE}"
done
echo 'Copying static files'
cp -rf 'docs/static/.' '../documentation_tmp'
REV=`git rev-parse HEAD`
DATE=`date +%Y-%m-%d_%H-%M-%S`
sed -i "s/REVISION/Build time: $DATE, Revision: $REV/g;" '../documentation_tmp/0_introduction.md'

cd ..
mkdir -p 'documentation'
if [ "$1" == 'html' ]; then
    echo 'Compiling to HTML'
    markdown-folder-to-html 'documentation_tmp'
    cp -rf '_documentation_tmp/.' 'documentation'
else
    echo 'Compiling to PDF'
    cd 'documentation_tmp'
    pandoc -r markdown_strict -o '../documentation/documentation.pdf' --filter '../automation/pandoc-svg.py' --toc --pdf-engine=xelatex -V geometry:margin=2cm -V documentclass=report `find . -name '*.md'`
    cd ..
fi

echo 'Cleaning up'
rm -rf 'documentation_tmp' '_documentation_tmp'
echo 'Documentation job ended'
