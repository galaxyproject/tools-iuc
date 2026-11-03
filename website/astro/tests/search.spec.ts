import { expect, test } from '@playwright/test';

test('generated Pagefind UI returns tool links', async ({ page }) => {
  await page.goto('tools/');
  const search = page.locator('#pagefind-search input');
  await expect(search).toBeVisible();
  await search.fill('fastqc');
  await expect(page.locator('#pagefind-search .pagefind-ui__result-link').first()).toHaveAttribute(
    'href', /\/tools\/[^/]+\/[^/]+\/[^/]+\//,
  );
});

test('search failure shows a usable fallback', async ({ page }) => {
  await page.route('**/pagefind/pagefind-ui.js', route => route.abort());
  await page.goto('tools/');
  await expect(page.locator('#pagefind-search')).toContainText('Full-text search is unavailable');
  await page.locator('#filter').fill('fastqc');
  await expect(page.locator('#tool-list > li:visible').first()).toContainText(/fastqc/i);
});

test('contributor profile separates Git contributions from creator credits', async ({ page }) => {
  await page.goto('contributors/mvdbeek/');
  const contributed = page.getByRole('heading', { name: 'Tools contributed to', exact: true }).locator('..');
  const authored = page.getByRole('heading', { name: 'Authored tools', exact: true }).locator('..');
  await expect(contributed.locator('a').first()).toHaveAttribute('href', /\/tools\//);
  expect(await contributed.locator('li').count()).toBeGreaterThan(await authored.locator('li').count());
  await expect(contributed.getByText('Example change').first()).toHaveAttribute('href', /github.com\/galaxyproject\/tools-iuc\/commit\//);
});
