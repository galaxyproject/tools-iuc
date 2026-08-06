#!/usr/bin/env node

/*
Use linkinator to crawl all internal links on the built Astro site.
Adapted from galaxy-hub's check-links.mjs.

Usage:
  npm run build && npm run links:internal
  or:
  npx astro preview & (sleep 5 && npm run links:internal)
*/
import { LinkChecker } from 'linkinator';
import fs from 'fs';

const port = process.env.PORT || 4321;
const base = process.env.BASE_PATH || '/tools-iuc';
const siteURL = `http://localhost:${port}${base}/`;
const outputFile = './broken-links.md';

const outTemplate = (pages, total, broken) => `
### Link summary of ${pages} pages checked

| Checked | Successful | Errors |
| ------: | ---------: | -----: |
|  ${total}  |    ${total - broken}    |   ${broken}   |

### Individual Details

`;

async function checkLinks() {
  console.log(`Checking internal links starting from ${siteURL}`);

  const brokenLinksByPage = {};
  let brokenCount = 0;
  let totalCount = 0;
  const pagesChecked = new Set();

  try {
    const checker = new LinkChecker();

    checker.on('link', (result) => {
      totalCount++;
      if (result.state === 'BROKEN') {
        brokenCount++;
        const parentPage = result.parent || siteURL;
        pagesChecked.add(parentPage);
        if (!brokenLinksByPage[parentPage]) {
          brokenLinksByPage[parentPage] = [];
        }
        brokenLinksByPage[parentPage].push({
          url: result.url,
          status: result.status,
          statusText: result.statusText,
        });
      }
    });

    checker.on('pagestart', (url) => {
      pagesChecked.add(url);
    });

    // Skip external links — we only check internal page links
    const skipLinkChecker = async (url) => {
      if (url.startsWith('http') && !url.startsWith(`http://localhost:${port}`)) {
        return true;
      }
      return false;
    };

    const result = await checker.check({
      path: siteURL,
      recurse: true,
      excludeExternalLinks: true,
      linksToSkip: skipLinkChecker,
      timeout: 30000,
      concurrency: 100,
    });

    console.log(`Completed checking ${result.links.length} links`);

    let markdownReport = '';
    for (const [page, links] of Object.entries(brokenLinksByPage)) {
      markdownReport += `#### ${page}\n`;
      for (const link of links) {
        markdownReport += `- [ ] ${link.url} (${link.status}: ${link.statusText || 'Unknown error'})\n`;
      }
      markdownReport += '\n';
    }

    const output = outTemplate(pagesChecked.size, totalCount, brokenCount) + markdownReport;
    fs.writeFileSync(outputFile, output);
    console.log(output);

    if (brokenCount > 0) {
      process.exitCode = 1;
    }
  } catch (error) {
    console.error('Error during link checking:', error);
    process.exitCode = 1;
  }
}

checkLinks();
